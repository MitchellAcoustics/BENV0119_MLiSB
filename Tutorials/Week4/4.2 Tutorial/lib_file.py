import numpy as np
import pandas as pd
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.neighbors         import KNeighborsRegressor  # The KNN regression model
from sklearn.model_selection   import cross_val_score      # The cross-validation function
from sklearn.preprocessing   import scale


def plot_cv_indices(cv, X, y, ax, n_splits, lw=10, cmap_cv = 'coolwarm'):
    """Create a sample plot for indices of a cross-validation object."""
    """Adapted from: https://scikit-learn.org/stable/auto_examples/model_selection/plot_cv_indices.html"""

    # Get colormap
    cmap_cv = plt.cm.get_cmap( cmap_cv )
    
    # Generate the training/testing visualizations for each CV split
    for ii, (tr, tt) in enumerate(cv.split(X=X, y=y)):
        # Fill in indices with the training/test groups
        indices = np.array([np.nan] * len(X))
        indices[tt] = 1
        indices[tr] = 0

        # Visualize the results
        ax.scatter(range(len(indices)), [ii + .5] * len(indices),
                   c=indices, marker='_', lw=lw, cmap=cmap_cv,
                   vmin=-.2, vmax=1.2)

    # Formatting
    yticklabels = list(range(n_splits))
    ax.set(yticks=np.arange(n_splits) + .5, yticklabels=yticklabels,
           xlabel='Sample index', ylabel="CV iteration",
           ylim=[n_splits+.2, -.2], xlim=[0, len(X)])
    ax.set_title('{}'.format(type(cv).__name__), fontsize=15)
    ax.legend([Patch(color=cmap_cv(.8)), Patch(color=cmap_cv(.02))],
              ['Validation set', 'Training set'], loc=(1.02, .8))
    return ax

def format_df(df, cmap, decimals = 3, vmin = None, vmax = None):
    if vmin is None: vmin = df.min().min()
    if vmax is None: vmax = df.max().max()
    return df.style.background_gradient(cmap=cmap, vmin = vmin, vmax = vmax).format('{:.%d}' %decimals) # Formatting

def map_radiation_class(Gt):
    cls = pd.Series(index = Gt.index, dtype = 'object')
    cls[:] = 'Very good'
    cls[Gt < 1400] = 'Good'
    cls[Gt < 1200] = 'Average'
    cls[Gt < 1000] = 'Bad'
    cls[Gt < 800 ] = 'Very bad'
    return cls

def cross_val_score_default_knn():
    t = 'tilted_radiation'
    f = ['roof_tilt', 'roof_aspect', 'roof_area', 'roof_x', 'roof_y', 'SIS',
         'SISDIR', 'DHI', 'horizon_S', 'horizon_SSW', 'horizon_SWW', 'horizon_W',
         'horizon_NWW', 'horizon_SSE', 'horizon_SEE', 'horizon_E', 'horizon_NEE',]

    d = pd.read_csv( 'training_sample.csv' , index_col = 0 )
    
    knn = KNeighborsRegressor() # Initiate model
    return cross_val_score(knn, scale(d[f]), d[t], cv = 5)