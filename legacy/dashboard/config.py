"""
MODFLOW Dashboard Configuration (Legacy)
========================================

NOTE: This is legacy code preserved for reference.

Professional configuration settings for the MODFLOW dashboard submodule.
"""

# =============================================================================
# VISUAL THEME CONFIGURATION
# =============================================================================

THEMES = {
    'professional': {
        'primary_color': '#1f4e79',
        'secondary_color': '#2e8b57',
        'accent_color': '#ff6b35',
        'background_color': '#f8f9fa',
        'text_color': '#2c3e50',
        'success_color': '#28a745',
        'warning_color': '#ffc107',
        'error_color': '#dc3545'
    },

    'academic': {
        'primary_color': '#2c3e50',
        'secondary_color': '#3498db',
        'accent_color': '#e74c3c',
        'background_color': '#ecf0f1',
        'text_color': '#2c3e50',
        'success_color': '#27ae60',
        'warning_color': '#f39c12',
        'error_color': '#e74c3c'
    },

    'colorblind_friendly': {
        'primary_color': '#0173b2',
        'secondary_color': '#029e73',
        'accent_color': '#cc78bc',
        'background_color': '#f8f9fa',
        'text_color': '#2c3e50',
        'success_color': '#029e73',
        'warning_color': '#fbafe4',
        'error_color': '#d55e00'
    },

    'high_contrast': {
        'primary_color': '#000000',
        'secondary_color': '#ffffff',
        'accent_color': '#ffff00',
        'background_color': '#ffffff',
        'text_color': '#000000',
        'success_color': '#008000',
        'warning_color': '#ffa500',
        'error_color': '#ff0000'
    }
}

# =============================================================================
# ANALYSIS CONFIGURATION
# =============================================================================

ANALYSIS_DEFAULTS = {
    'time_series': {
        'default_variables': ['Qaquifer', 'stage'],
        'default_reaches': [0, 39, 78],  # First, middle, last
        'smoothing_window': 15,
        'gaussian_sigma': 1.5,
        'spectral_nperseg': 256
    },

    'river_analysis': {
        'flow_components': ['Qin', 'Qout', 'Qaquifer'],
        'hydraulic_variables': ['stage', 'velocity', 'hydraulic_radius'],
        'exchange_threshold': 1e-6,  # Minimum exchange flow to consider
        'representative_reaches': [0, 20, 40, 60, 78]
    },

    'statistics': {
        'confidence_level': 0.95,
        'outlier_method': 'IQR',  # 'IQR' or 'zscore'
        'outlier_threshold': 1.5,  # For IQR method
        'correlation_method': 'pearson',  # 'pearson', 'spearman', 'kendall'
    },

    'export': {
        'default_dpi': 300,
        'figure_formats': ['PNG', 'PDF', 'SVG', 'EPS'],
        'data_formats': ['CSV', 'Excel', 'JSON', 'HDF5'],
        'compression_formats': ['ZIP', 'TAR.GZ', '7Z']
    }
}

# =============================================================================
# PLOTLY CONFIGURATION
# =============================================================================

PLOTLY_CONFIG = {
    'displayModeBar': True,
    'displaylogo': False,
    'modeBarButtonsToRemove': [
        'pan2d', 'select2d', 'lasso2d', 'autoScale2d'
    ],
    'toImageButtonOptions': {
        'format': 'png',
        'filename': 'modflow_analysis',
        'height': 800,
        'width': 1200,
        'scale': 2
    }
}

# =============================================================================
# DASHBOARD LAYOUT CONFIGURATION
# =============================================================================

LAYOUT_CONFIG = {
    'sidebar_width': 300,
    'main_content_padding': '2rem',
    'header_height': 120,
    'footer_height': 60,
    'card_border_radius': '8px',
    'card_shadow': '0 2px 10px rgba(0,0,0,0.1)'
}

# =============================================================================
# PERFORMANCE CONFIGURATION
# =============================================================================

PERFORMANCE_CONFIG = {
    'cache_ttl': 3600,  # Cache time-to-live in seconds
    'max_data_points': 10000,  # Maximum points for interactive plots
    'animation_fps': 10,  # Default animation frame rate
    'chunk_size': 1000,  # Data processing chunk size
    'memory_limit_mb': 512  # Memory limit for large datasets
}

# =============================================================================
# SCIENTIFIC CONFIGURATION
# =============================================================================

SCIENTIFIC_CONFIG = {
    'units': {
        'length': 'm',
        'time': 's',
        'discharge': 'm3/s',
        'head': 'm',
        'conductivity': 'm/s',
        'area': 'm2',
        'volume': 'm3'
    },

    'conversion_factors': {
        'time': {
            'seconds_to_minutes': 1/60,
            'seconds_to_hours': 1/3600,
            'seconds_to_days': 1/86400
        },
        'discharge': {
            'm3_s_to_m3_day': 86400,
            'm3_s_to_l_s': 1000,
            'm3_s_to_cfs': 35.3147
        }
    },

    'precision': {
        'coordinates': 2,
        'flows': 4,
        'heads': 3,
        'statistics': 4,
        'percentages': 2
    }
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_theme(theme_name='professional'):
    """Get theme configuration by name."""
    return THEMES.get(theme_name, THEMES['professional'])

def get_analysis_defaults(analysis_type):
    """Get default parameters for analysis type."""
    return ANALYSIS_DEFAULTS.get(analysis_type, {})

def format_number(value, precision_type='default'):
    """Format number according to scientific configuration."""
    precision = SCIENTIFIC_CONFIG['precision'].get(precision_type, 4)
    return f"{value:.{precision}f}"
