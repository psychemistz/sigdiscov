#!/usr/bin/env python
"""
genomewide_interaction_v9_auto_optimized.py - AUTO-OPTIMIZED FACTOR-SPECIFIC VERSION

Intelligent batch size selection based on number of factors and available GPU memory
Automatically optimizes for best performance based on problem size

Key features:
1. Auto-detects optimal batch size based on number of factors
2. Reads factor list from TSV file (optional)
3. Fixes GPU array conversion issues
4. Maximizes GPU utilization for any factor count
"""

import os
import argparse
from typing import List, Tuple, Dict, Optional, Any
import numpy as np
import pandas as pd
from pathlib import Path
import time
from datetime import datetime, timedelta
import anndata as ad
from scipy import sparse
import h5py
from tqdm import tqdm
import gc
import warnings
import logging
from dataclasses import dataclass, field
import psutil
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import threading
from queue import Queue
import numba as nb
import gzip

warnings.filterwarnings('ignore')

# ============================================================================
# OPTIMIZED CONFIGURATION WITH AUTO BATCH SIZING
# ============================================================================

@dataclass
class AutoOptimizedConfig:
    """Configuration with intelligent batch size selection"""
    
    # Core parameters
    RANDOM_SEED: int = 42
    RADII: List[float] = field(default_factory=lambda: [10, 20, 30, 50, 100, 200, 300, 500])
    MIN_EXPR_QUANTILE: float = 0.25
    MIN_CELLS_THRESHOLD: int = 10
    MIN_CONNECTIONS_THRESHOLD: int = 10
    
    # Optimization parameters (will be auto-adjusted)
    USE_FLOAT16: bool = True
    FACTOR_BATCH_SIZE: int = None  # Will be set automatically
    TARGET_BATCH_SIZE: int = 1000
    MAX_CHUNK_SIZE: int = 5000
    
    # Annular parameters
    USE_ANNULAR_EDGES: bool = False
    ANNULAR_WIDTH: float = 20.0
    
    # Output parameters
    USE_HDF5: bool = True
    COMPRESSION: str = 'gzip'
    COMPRESSION_LEVEL: int = 1
    
    # Memory management
    CLEAR_GPU_EVERY_N: int = 10
    USE_MEMORY_MAPPING: bool = True
    
    # Auto-optimization parameters
    GPU_MEMORY_FRACTION: float = 0.8  # Use 80% of available GPU memory
    MIN_BATCH_SIZE: int = 50  # Minimum batch size
    MAX_BATCH_SIZE: int = 2000  # Maximum batch size
    
    # Parallel processing
    N_CPU_WORKERS: int = 4
    
    # GPU_AVAILABLE_MEM will be set dynamically based on GPU type
    GPU_AVAILABLE_MEM: float = 30e9  # Default, will be updated
    
    def get_sigma(self, radius: float) -> float:
        return radius / 3.0
    
    def get_inner_radius(self, outer_radius: float) -> Optional[float]:
        if not self.USE_ANNULAR_EDGES:
            return None
        return max(0, outer_radius - self.ANNULAR_WIDTH)
    
    @classmethod
    def generate_regular_radii(cls, max_distance: float = 5000, spacing: float = 20) -> List[float]:
        return list(range(int(spacing), int(max_distance) + 1, int(spacing)))
    
    def calculate_optimal_batch_size(self, n_factors: int, n_cells: int, n_targets: int, 
                                    available_gpu_memory: float, logger: logging.Logger) -> int:
        """
        Calculate optimal batch size based on problem dimensions and GPU memory
        """
        # Estimate memory per factor in batch (in GB)
        # Memory needed: factor data + masks + intermediate computations
        bytes_per_element = 2 if self.USE_FLOAT16 else 4
        
        # Memory components per factor:
        # - Factor expression: n_cells * bytes_per_element
        # - Quantile mask: n_cells bytes (boolean)
        # - Spatial lags subset: n_cells * n_targets * bytes_per_element (worst case)
        # - Intermediate computations: ~2x buffer
        
        memory_per_factor = (
            n_cells * bytes_per_element +  # Factor data
            n_cells +  # Mask
            n_cells * bytes_per_element * 10  # Intermediate (conservative estimate)
        ) / 1e9  # Convert to GB
        
        # Available memory (in GB)
        usable_memory = available_gpu_memory * self.GPU_MEMORY_FRACTION
        
        # Calculate max factors that fit in memory
        max_factors_in_memory = int(usable_memory / memory_per_factor)
        
        # Determine optimal batch size
        if n_factors <= self.MAX_BATCH_SIZE and n_factors <= max_factors_in_memory:
            # Process all factors at once if they fit
            optimal_batch = n_factors
            logger.info(f"✓ All {n_factors} factors fit in memory - single batch processing")
        else:
            # Use the smaller of max_memory or max_batch_size
            optimal_batch = min(max_factors_in_memory, self.MAX_BATCH_SIZE)
            optimal_batch = max(optimal_batch, self.MIN_BATCH_SIZE)  # Ensure minimum
            n_batches = (n_factors + optimal_batch - 1) // optimal_batch
            logger.info(f"Will process {n_factors} factors in {n_batches} batches of {optimal_batch}")
        
        logger.info(f"Memory estimate: {memory_per_factor:.3f}GB per factor, "
                   f"{usable_memory:.1f}GB available")
        
        return optimal_batch

# ============================================================================
# NUMBA-ACCELERATED FUNCTIONS
# ============================================================================

@nb.njit(parallel=True, fastmath=True)
def compute_distances_chunked(sender_coords: np.ndarray, receiver_coords: np.ndarray,
                              radius: float, inner_radius: float = 0.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Numba-accelerated distance computation"""
    n_senders = sender_coords.shape[0]
    n_receivers = receiver_coords.shape[0]
    
    rows = []
    cols = []
    dists = []
    
    radius_sq = radius * radius
    inner_radius_sq = inner_radius * inner_radius if inner_radius > 0 else 0.0
    
    for i in nb.prange(n_senders):
        for j in range(n_receivers):
            dx = sender_coords[i, 0] - receiver_coords[j, 0]
            dy = sender_coords[i, 1] - receiver_coords[j, 1]
            dist_sq = dx * dx + dy * dy
            
            if dist_sq <= radius_sq and dist_sq > inner_radius_sq:
                rows.append(i)
                cols.append(j)
                dists.append(np.sqrt(dist_sq))
    
    return np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32), np.array(dists, dtype=np.float32)

# ============================================================================
# FACTOR LIST READER
# ============================================================================

def read_factor_list(factor_file_path: str, logger: logging.Logger) -> Optional[List[str]]:
    """Read factor gene names from TSV file"""
    if not factor_file_path:
        return None
        
    logger.info(f"Reading factor list from {factor_file_path}")
    
    try:
        if factor_file_path.endswith('.gz'):
            with gzip.open(factor_file_path, 'rt') as f:
                df = pd.read_csv(f, sep='\t', nrows=1)
        else:
            df = pd.read_csv(factor_file_path, sep='\t', nrows=1)
        
        factor_names = df.columns[1:].tolist()
        logger.info(f"Found {len(factor_names)} factors in file")
        logger.info(f"First 5 factors: {factor_names[:5]}")
        
        return factor_names
        
    except Exception as e:
        logger.error(f"Error reading factor file: {e}")
        raise

# ============================================================================
# OPTIMIZED WEIGHT MATRIX BUILDER
# ============================================================================

class OptimizedWeightBuilder:
    """Optimized weight matrix construction"""
    
    def __init__(self, cp, config: AutoOptimizedConfig):
        self.cp = cp
        self.config = config
        self.cache = {}
        
    def build_weight_matrix_ultra_fast(self, sender_coords_gpu, receiver_coords_gpu,
                                      radius: float, inner_radius: Optional[float] = None):
        """Ultra-fast weight matrix construction"""
        
        n_senders = len(sender_coords_gpu)
        n_receivers = len(receiver_coords_gpu)
        
        dtype = self.cp.float16 if self.config.USE_FLOAT16 else self.cp.float32
        
        mem_per_chunk = n_receivers * 4 * (2 if self.config.USE_FLOAT16 else 4)
        # Dynamically set available memory based on GPU
        # Will be properly set when GPU is initialized
        available_mem = getattr(self.config, 'GPU_AVAILABLE_MEM', 30e9)
        chunk_size = min(n_senders, int(available_mem / mem_per_chunk / 4))
        
        sigma = self.config.get_sigma(radius)
        gaussian_factor = -1.0 / (2 * sigma * sigma)
        radius_sq = radius * radius
        inner_radius_sq = (inner_radius * inner_radius) if inner_radius else 0
        
        rows_list = []
        cols_list = []
        weights_list = []
        
        for chunk_start in range(0, n_senders, chunk_size):
            chunk_end = min(chunk_start + chunk_size, n_senders)
            
            sender_chunk = sender_coords_gpu[chunk_start:chunk_end].astype(dtype)
            
            diff_x = sender_chunk[:, 0:1] - receiver_coords_gpu[:, 0].astype(dtype).T
            diff_y = sender_chunk[:, 1:2] - receiver_coords_gpu[:, 1].astype(dtype).T
            
            dist_sq = diff_x * diff_x + diff_y * diff_y
            
            mask = dist_sq <= radius_sq
            if inner_radius:
                mask = mask & (dist_sq > inner_radius_sq)
            
            weights_chunk = self.cp.exp(dist_sq * gaussian_factor) * mask
            
            nonzero_mask = weights_chunk > 1e-6
            if self.cp.any(nonzero_mask):
                rows, cols = self.cp.where(nonzero_mask)
                rows_list.append(rows + chunk_start)
                cols_list.append(cols)
                weights_list.append(weights_chunk[nonzero_mask].astype(self.cp.float32))
            
            del diff_x, diff_y, dist_sq, mask, weights_chunk
            
        if rows_list:
            all_rows = self.cp.concatenate(rows_list)
            all_cols = self.cp.concatenate(cols_list)
            all_weights = self.cp.concatenate(weights_list)
            
            W_sparse = self.cp.sparse.coo_matrix(
                (all_weights, (all_rows, all_cols)),
                shape=(n_senders, n_receivers),
                dtype=self.cp.float32
            ).tocsr()
            
            row_sums = self.cp.array(W_sparse.sum(axis=1)).ravel()
            row_sums_inv = self.cp.where(row_sums > 0, 1.0 / row_sums, 0)
            D_inv = self.cp.sparse.diags(row_sums_inv, dtype=self.cp.float32, format='csr')
            W_normalized = D_inv @ W_sparse
            
            return W_normalized, float(W_sparse.sum())
        else:
            return None, 0.0

# ============================================================================
# OPTIMIZED MORAN'S I COMPUTATION
# ============================================================================

class MoransIMatrixComputer:
    """Optimized Moran's I computation with dynamic batching"""
    
    def __init__(self, cp, config: AutoOptimizedConfig):
        self.cp = cp
        self.config = config
        
    def compute_morans_matrix_batched(self, 
                                     all_factor_expr_gpu,
                                     all_target_expr_gpu,
                                     W_normalized,
                                     quantile_masks_gpu):
        """
        Compute Moran's I matrix with intelligent batching
        """
        
        n_factors = all_factor_expr_gpu.shape[1]
        n_targets = all_target_expr_gpu.shape[1]
        
        # Pre-compute spatial lags for all targets
        spatial_lags_all = W_normalized @ all_target_expr_gpu
        
        dtype = self.cp.float16 if self.config.USE_FLOAT16 else self.cp.float32
        morans_matrix = self.cp.zeros((n_factors, n_targets), dtype=dtype)
        
        # Use configured batch size (which was auto-optimized)
        factor_batch_size = self.config.FACTOR_BATCH_SIZE
        
        # If processing all at once, use optimized path
        if factor_batch_size >= n_factors:
            # Single batch - process all factors
            for i in range(n_factors):
                high_expr_mask = quantile_masks_gpu[i]
                
                if self.cp.sum(high_expr_mask) < self.config.MIN_CELLS_THRESHOLD:
                    continue
                
                high_expr_indices = self.cp.where(high_expr_mask)[0]
                factor_high = all_factor_expr_gpu[high_expr_indices, i]
                spatial_lags_high = spatial_lags_all[high_expr_indices, :]
                
                factor_norm = self.cp.linalg.norm(factor_high)
                if factor_norm > 0:
                    factor_normalized = factor_high / factor_norm
                    
                    correlations = self.cp.dot(factor_normalized, spatial_lags_high)
                    spatial_norms = self.cp.linalg.norm(spatial_lags_high, axis=0)
                    valid_mask = spatial_norms > 0
                    
                    correlations[valid_mask] /= spatial_norms[valid_mask]
                    morans_matrix[i, :] = correlations
        else:
            # Multiple batches needed
            for f_start in range(0, n_factors, factor_batch_size):
                f_end = min(f_start + factor_batch_size, n_factors)
                batch_size = f_end - f_start
                
                factor_batch = all_factor_expr_gpu[:, f_start:f_end]
                masks_batch = quantile_masks_gpu[f_start:f_end, :]
                
                for i in range(batch_size):
                    factor_idx = f_start + i
                    high_expr_mask = masks_batch[i]
                    
                    if self.cp.sum(high_expr_mask) < self.config.MIN_CELLS_THRESHOLD:
                        continue
                    
                    high_expr_indices = self.cp.where(high_expr_mask)[0]
                    factor_high = factor_batch[high_expr_indices, i]
                    spatial_lags_high = spatial_lags_all[high_expr_indices, :]
                    
                    factor_norm = self.cp.linalg.norm(factor_high)
                    if factor_norm > 0:
                        factor_normalized = factor_high / factor_norm
                        
                        correlations = self.cp.dot(factor_normalized, spatial_lags_high)
                        spatial_norms = self.cp.linalg.norm(spatial_lags_high, axis=0)
                        valid_mask = spatial_norms > 0
                        
                        correlations[valid_mask] /= spatial_norms[valid_mask]
                        morans_matrix[factor_idx, :] = correlations
        
        return morans_matrix

# ============================================================================
# HDF5 OUTPUT HANDLER
# ============================================================================

class HDF5OutputHandler:
    """HDF5 output handler"""
    
    def __init__(self, output_path: Path, config: AutoOptimizedConfig):
        self.output_path = output_path
        self.config = config
        self.file = None
        self.lock = threading.Lock()
        
    def __enter__(self):
        self.file = h5py.File(self.output_path, 'w')
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.file:
            self.file.close()
            
    def initialize_dataset(self, n_pairs: int, n_radii: int, n_factors: int, n_targets: int,
                          cell_types: List[str], radii: List[float], 
                          factor_names: List[str], target_names: List[str]):
        """Initialize HDF5 structure"""
        
        dtype = np.float16 if self.config.USE_FLOAT16 else np.float32
        
        self.file.create_dataset(
            'morans_i_matrices',
            shape=(n_pairs, n_radii, n_factors, n_targets),
            dtype=dtype,
            chunks=(1, 1, min(100, n_factors), min(100, n_targets)),
            compression=self.config.COMPRESSION,
            compression_opts=self.config.COMPRESSION_LEVEL
        )
        
        self.file.attrs['n_pairs'] = n_pairs
        self.file.attrs['n_radii'] = n_radii
        self.file.attrs['n_factors'] = n_factors
        self.file.attrs['n_targets'] = n_targets
        self.file.attrs['factor_batch_size'] = self.config.FACTOR_BATCH_SIZE
        
        self.file.create_dataset('cell_types', data=np.array(cell_types, dtype='S'))
        self.file.create_dataset('radii', data=radii)
        self.file.create_dataset('factor_genes', data=np.array(factor_names, dtype='S'))
        self.file.create_dataset('target_genes', data=np.array(target_names, dtype='S'))
        
        pairs = []
        for i, sender in enumerate(cell_types):
            for j, receiver in enumerate(cell_types):
                if i != j:
                    pairs.append(f"{sender}->{receiver}")
        self.file.create_dataset('pairs', data=np.array(pairs, dtype='S'))
        
    def write_matrix(self, pair_idx: int, radius_idx: int, morans_matrix: np.ndarray):
        """Write matrix to HDF5"""
        with self.lock:
            self.file['morans_i_matrices'][pair_idx, radius_idx, :, :] = morans_matrix
            self.file.flush()

# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

def run_matrix_analysis(adata: ad.AnnData,
                       output_dir: Path,
                       config: AutoOptimizedConfig,
                       cell_type_col: str,
                       gpu_id: str,
                       logger: logging.Logger,
                       factor_names: Optional[List[str]] = None,
                       gpu_type: str = 'auto') -> Dict[str, Any]:
    """
    Run matrix analysis with auto-optimized batch sizing
    
    Args:
        gpu_type: 'auto', 'a100', 'v100', or memory size in GB (e.g., '80')
    """
    
    # Setup GPU
    os.environ['CUDA_VISIBLE_DEVICES'] = gpu_id
    try:
        import cupy as cp
        mempool = cp.get_default_memory_pool()
        
        # Determine GPU memory based on type or auto-detect
        if gpu_type == 'auto':
            try:
                # Try to auto-detect GPU memory
                mem_info = cp.cuda.runtime.memGetInfo()
                total_memory_bytes = mem_info[1]
                gpu_memory_gb = int(total_memory_bytes / (1024**3) * 0.9)  # Use 90% of detected
                logger.info(f"Auto-detected GPU memory: {total_memory_bytes/(1024**3):.1f}GB")
            except:
                # Default to conservative estimate
                gpu_memory_gb = 35
                logger.warning("Could not auto-detect GPU memory, using 35GB default")
        elif gpu_type.lower() == 'a100':
            gpu_memory_gb = 75  # Use 75GB out of 80GB
            logger.info("GPU type: A100 (80GB total)")
        elif gpu_type.lower() == 'v100':
            gpu_memory_gb = 35  # Use 35GB out of 40GB
            logger.info("GPU type: V100 (40GB total)")
        else:
            # Assume it's a number representing GB
            try:
                total_gb = float(gpu_type)
                gpu_memory_gb = int(total_gb * 0.9)  # Use 90% of specified
                logger.info(f"GPU memory specified: {total_gb}GB total")
            except:
                gpu_memory_gb = 35
                logger.warning(f"Could not parse GPU type '{gpu_type}', using 35GB default")
        
        # Set memory pool limit
        mempool.set_limit(size=gpu_memory_gb * 1024**3)
        logger.info(f"✓ GPU initialized with {gpu_memory_gb}GB memory limit")
        
        # Store available memory in config for weight builder
        config.GPU_AVAILABLE_MEM = gpu_memory_gb * 0.9 * 1e9  # Use 90% for computations
    except Exception as e:
        logger.error(f"GPU setup failed: {e}")
        return {'success': False}
    
    # Filter genes
    logger.info("Filtering genes...")
    n_cells = adata.shape[0]
    if sparse.issparse(adata.X):
        X_csc = adata.X.tocsc()
        gene_counts = np.diff(X_csc.indptr)
    else:
        gene_counts = (adata.X > 0).sum(axis=0)
    
    min_cells = max(10, n_cells * 0.001)
    keep_genes = gene_counts >= min_cells
    adata_filtered = adata[:, keep_genes].copy()
    
    n_genes = adata_filtered.shape[1]
    logger.info(f"Genes: {n_genes}/{adata.shape[1]} kept (min {min_cells:.0f} cells)")
    
    # Determine factors
    gene_names = adata_filtered.var_names.tolist()
    
    if factor_names is not None:
        logger.info("Running FACTOR-SPECIFIC analysis")
        factor_indices = []
        factor_names_found = []
        
        for factor in factor_names:
            if factor in gene_names:
                idx = gene_names.index(factor)
                factor_indices.append(idx)
                factor_names_found.append(factor)
        
        n_factors = len(factor_indices)
        logger.info(f"Found {n_factors}/{len(factor_names)} factors in gene list")
        
        if n_factors == 0:
            logger.error("No factors found in gene list!")
            return {'success': False}
    else:
        logger.info("Running FULL analysis (all genes as factors)")
        factor_indices = list(range(n_genes))
        factor_names_found = gene_names
        n_factors = n_genes
    
    # AUTO-OPTIMIZE BATCH SIZE
    optimal_batch_size = config.calculate_optimal_batch_size(
        n_factors=n_factors,
        n_cells=n_cells,
        n_targets=n_genes,
        available_gpu_memory=gpu_memory_gb,
        logger=logger
    )
    config.FACTOR_BATCH_SIZE = optimal_batch_size
    logger.info(f"✓ Optimal batch size set to: {optimal_batch_size}")
    
    # Get cell types
    cell_types = []
    for ct in adata.obs[cell_type_col].unique():
        if (adata.obs[cell_type_col] == ct).sum() >= config.MIN_CELLS_THRESHOLD:
            cell_types.append(ct)
    
    # Generate pairs
    pairs = []
    for i, sender in enumerate(cell_types):
        for j, receiver in enumerate(cell_types):
            if i != j:
                pairs.append((sender, receiver))
    
    n_pairs = len(pairs)
    logger.info(f"Cell type pairs: {n_pairs} (from {len(cell_types)} types)")
    
    radii = config.RADII
    n_radii = len(radii)
    
    # Pre-compute normalized expression
    logger.info("Pre-computing normalized expression matrices...")
    
    if sparse.issparse(adata_filtered.X):
        logger.info("Converting sparse to dense...")
        X_dense = adata_filtered.X.toarray().astype(np.float32)
    else:
        X_dense = adata_filtered.X.astype(np.float32)
    
    X_mean = X_dense.mean(axis=0, keepdims=True)
    X_std = X_dense.std(axis=0, keepdims=True)
    X_normalized = (X_dense - X_mean) / (X_std + 1e-10)
    
    # Extract factor data
    if factor_names is not None:
        X_factors = X_dense[:, factor_indices]
        X_normalized_factors = X_normalized[:, factor_indices]
    else:
        X_factors = X_dense
        X_normalized_factors = X_normalized
    
    # Compute quantile thresholds
    logger.info(f"Computing {config.MIN_EXPR_QUANTILE:.0%} quantile thresholds for {n_factors} factors...")
    quantile_thresholds_factors = np.quantile(X_factors, config.MIN_EXPR_QUANTILE, axis=0)
    
    # Convert to GPU
    quantile_thresholds_factors_gpu = cp.asarray(quantile_thresholds_factors, dtype=cp.float32)
    logger.info("Quantile thresholds moved to GPU")
    
    # Cache cell data
    logger.info("Caching cell coordinates...")
    cell_coords_cache = {}
    cell_expr_cache = {}
    cell_factor_cache = {}
    
    for ct in cell_types:
        mask = (adata_filtered.obs[cell_type_col] == ct).values
        coords = adata_filtered.obsm['spatial'][mask].astype(np.float32)
        cell_coords_cache[ct] = {
            'coords': coords,
            'mask': mask,
            'n_cells': mask.sum()
        }
        cell_expr_cache[ct] = X_normalized[mask, :]
        cell_factor_cache[ct] = X_normalized_factors[mask, :]
    
    # Initialize output
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    if factor_names is not None:
        output_file = output_dir / f'morans_matrix_factors_{timestamp}.h5'
    else:
        output_file = output_dir / f'morans_matrix_{timestamp}.h5'
    
    # Initialize helpers
    weight_builder = OptimizedWeightBuilder(cp, config)
    morans_computer = MoransIMatrixComputer(cp, config)
    
    # Main processing
    logger.info("="*80)
    if factor_names is not None:
        logger.info(f"Starting FACTOR-SPECIFIC computation:")
        logger.info(f"  Factors: {n_factors} (vs {n_genes} total genes)")
        logger.info(f"  Batch size: {config.FACTOR_BATCH_SIZE}")
        logger.info(f"  Batches needed: {(n_factors + config.FACTOR_BATCH_SIZE - 1) // config.FACTOR_BATCH_SIZE}")
        logger.info(f"  Reduction: {n_genes/n_factors:.1f}x fewer computations")
    else:
        logger.info(f"Starting FULL matrix computation:")
        logger.info(f"  Batch size: {config.FACTOR_BATCH_SIZE}")
    logger.info(f"  Pairs: {n_pairs}")
    logger.info(f"  Radii: {n_radii}")
    logger.info(f"  Matrix size: {n_pairs} × {n_radii} × {n_factors} × {n_genes}")
    logger.info(f"  Total values: {n_pairs * n_radii * n_factors * n_genes:,}")
    logger.info("="*80)
    
    start_time = time.time()
    
    with HDF5OutputHandler(output_file, config) as h5_handler:
        h5_handler.initialize_dataset(
            n_pairs, n_radii, n_factors, n_genes,
            cell_types, radii,
            factor_names_found,
            adata_filtered.var_names.tolist()
        )
        
        total_operations = n_pairs * n_radii
        with tqdm(total=total_operations, desc="Computing matrices") as pbar:
            
            for radius_idx, radius in enumerate(radii):
                inner_radius = config.get_inner_radius(radius) if config.USE_ANNULAR_EDGES else None
                
                for pair_idx, (sender, receiver) in enumerate(pairs):
                    
                    sender_info = cell_coords_cache[sender]
                    receiver_info = cell_coords_cache[receiver]
                    
                    sender_coords = sender_info['coords']
                    receiver_coords = receiver_info['coords']
                    sender_factors = cell_factor_cache[sender]
                    receiver_expr = cell_expr_cache[receiver]
                    
                    sender_coords_gpu = cp.asarray(sender_coords)
                    receiver_coords_gpu = cp.asarray(receiver_coords)
                    
                    W_normalized, S0 = weight_builder.build_weight_matrix_ultra_fast(
                        sender_coords_gpu, receiver_coords_gpu,
                        radius, inner_radius
                    )
                    
                    if W_normalized is None or W_normalized.nnz < config.MIN_CONNECTIONS_THRESHOLD:
                        zeros = np.zeros((n_factors, n_genes), dtype=np.float16)
                        h5_handler.write_matrix(pair_idx, radius_idx, zeros)
                        pbar.update(1)
                        continue
                    
                    sender_factors_gpu = cp.asarray(sender_factors, dtype=cp.float32)
                    receiver_expr_gpu = cp.asarray(receiver_expr, dtype=cp.float32)
                    
                    sender_quantile_masks = sender_factors_gpu > quantile_thresholds_factors_gpu[cp.newaxis, :]
                    sender_quantile_masks = sender_quantile_masks.T
                    
                    morans_matrix = morans_computer.compute_morans_matrix_batched(
                        sender_factors_gpu,
                        receiver_expr_gpu,
                        W_normalized,
                        sender_quantile_masks
                    )
                    
                    morans_matrix_cpu = cp.asnumpy(morans_matrix)
                    h5_handler.write_matrix(pair_idx, radius_idx, morans_matrix_cpu)
                    
                    del W_normalized, sender_factors_gpu, receiver_expr_gpu, sender_quantile_masks
                    
                    pbar.update(1)
                    
                    elapsed = time.time() - start_time
                    rate = (pair_idx + 1 + radius_idx * n_pairs) / elapsed
                    eta_seconds = (total_operations - (pair_idx + 1 + radius_idx * n_pairs)) / rate
                    eta_str = str(timedelta(seconds=int(eta_seconds)))
                    pbar.set_postfix_str(f"R={radius:.0f}μm | ETA={eta_str}")
                
                if (radius_idx + 1) % config.CLEAR_GPU_EVERY_N == 0:
                    mempool.free_all_blocks()
                    cp._default_memory_pool.free_all_blocks()
                    gc.collect()
    
    runtime = (time.time() - start_time) / 60
    
    logger.info("="*80)
    logger.info(f"✓ COMPLETE - Runtime: {runtime:.1f} minutes ({runtime/60:.1f} hours)")
    logger.info(f"Output: {output_file}")
    logger.info(f"Matrix shape: {n_pairs} × {n_radii} × {n_factors} × {n_genes}")
    logger.info(f"Total values: {n_pairs * n_radii * n_factors * n_genes:,}")
    logger.info(f"Speed: {(n_pairs * n_radii * n_factors * n_genes) / (runtime * 60):.0f} values/sec")
    
    if factor_names is not None:
        reduction = (n_genes * n_genes) / (n_factors * n_genes)
        logger.info(f"Computation reduction: {reduction:.1f}x vs full matrix")
    
    logger.info("="*80)
    
    return {
        'success': True,
        'output_file': output_file,
        'runtime_minutes': runtime,
        'shape': (n_pairs, n_radii, n_factors, n_genes)
    }

# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def setup_logging(output_dir: Path) -> logging.Logger:
    """Setup logging"""
    logger = logging.getLogger('spatial_morans_matrix')
    logger.setLevel(logging.DEBUG)
    
    if logger.handlers:
        logger.handlers.clear()
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    fh = logging.FileHandler(output_dir / f'analysis_log_{timestamp}.log')
    fh.setLevel(logging.DEBUG)
    
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    
    logger.addHandler(fh)
    logger.addHandler(ch)
    
    return logger

def main():
    parser = argparse.ArgumentParser(description='Auto-Optimized Genome-Wide Spatial Morans I Matrix')
    
    # All original arguments
    parser.add_argument('--data', type=str, required=True, help='Path to h5ad file')
    parser.add_argument('--output', type=str, required=True, help='Output directory')
    parser.add_argument('--cell-type-col', type=str, default='cell_type')
    parser.add_argument('--min-expr-quantile', type=float, default=0.25)
    parser.add_argument('--regular-radii', action='store_true')
    parser.add_argument('--max-distance', type=float, default=5000)
    parser.add_argument('--spacing', type=float, default=20)
    parser.add_argument('--radii', type=float, nargs='+', default=None)
    parser.add_argument('--use-annular-edges', action='store_true')
    parser.add_argument('--annular-width', type=float, default=20)
    parser.add_argument('--gpu', type=str, default='0')
    parser.add_argument('--use-float16', action='store_true', help='Use half precision')
    
    # GPU type specification
    parser.add_argument('--gpu-type', type=str, default='auto',
                       help='GPU type: auto, a100, v100, or memory in GB (default: auto-detect)')
    
    # Modified to support auto-optimization
    parser.add_argument('--factor-batch-size', type=int, default=None,
                       help='Override auto batch size (default: automatic based on factor count)')
    
    # Factor file for factor-specific computation
    parser.add_argument('--factors', type=str, default=None, 
                       help='Optional: Path to factor TSV file for factor-specific computation')
    
    args = parser.parse_args()
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    logger = setup_logging(output_dir)
    
    logger.info("="*80)
    logger.info("AUTO-OPTIMIZED SPATIAL MORAN'S I MATRIX")
    if args.gpu_type != 'auto':
        logger.info(f"GPU Type: {args.gpu_type.upper()}")
    logger.info("="*80)
    
    # Load data
    logger.info(f"Loading data from {args.data}...")
    adata = ad.read_h5ad(args.data)
    logger.info(f"Loaded: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    # Read factor list if provided
    factor_names = None
    if args.factors:
        factor_names = read_factor_list(args.factors, logger)
    
    # Configure
    config = AutoOptimizedConfig()
    config.MIN_EXPR_QUANTILE = args.min_expr_quantile
    config.USE_ANNULAR_EDGES = args.use_annular_edges
    config.ANNULAR_WIDTH = args.annular_width
    config.USE_FLOAT16 = args.use_float16
    
    # Override batch size if specified
    if args.factor_batch_size is not None:
        config.FACTOR_BATCH_SIZE = args.factor_batch_size
        logger.info(f"Using manual batch size: {args.factor_batch_size}")
    # Otherwise will be auto-calculated in run_matrix_analysis
    
    if args.regular_radii:
        config.RADII = config.generate_regular_radii(args.max_distance, args.spacing)
        logger.info(f"Radii: {len(config.RADII)} distances from {args.spacing} to {args.max_distance}")
    elif args.radii:
        config.RADII = args.radii
    
    # Run analysis
    result = run_matrix_analysis(
        adata, output_dir, config,
        args.cell_type_col, args.gpu, logger,
        factor_names=factor_names,
        gpu_type=args.gpu_type
    )
    
    if result['success']:
        logger.info(f"✓ Analysis completed successfully!")
        logger.info(f"Output saved to: {result['output_file']}")
        return 0
    else:
        logger.error("Analysis failed")
        return 1

if __name__ == "__main__":
    import sys
    sys.exit(main())
