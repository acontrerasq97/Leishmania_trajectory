import pandas as pd
import scanpy as sc

#rbps fropm haptomonads trajectory
rbps_1 = pd.read_csv('leishmania_trajectory/results/wot/predictive_rbps/RBPS_haptomonads_9_5.csv')
top_rbps_1 = rbps_1['gene_id']

#rbps fropm late-metacyclics trajectory
rbps_2 = pd.read_csv('leishmania_trajectory/results/wot/predictive_rbps/RBPS_late_metacyclics_9_5.csv')
top_rbps_2 = rbps_2['gene_id']

#export expression matrix 
adata = sc.read_h5ad('leishmania_trajectory/data/processed/ExprMatrix.preprocessed.h5ad')

#EM haptomonads
expr_haptomonads = adata[:, top_rbps_1].to_df()
expr_haptomonads['day'] = adata.obs['day'].values
expr_haptomonads.to_csv('leishmania_trajectory/data/processed/ExprMatrix_rbps_haptomonads.csv', index_label='index')

#EM haptomonads
expr_metacyclics = adata[:, top_rbps_2].to_df()
expr_metacyclics['day'] = adata.obs['day'].values
expr_metacyclics.to_csv('leishmania_trajectory/data/processed/ExprMatrix_rbps_late_metacyclics.csv', index_label='index')

#fates to dataframe
#fates = sc.read_h5ad('leishmania_trajectory/results/wot/fate_matrices/fate_matrices_model.h5ad')
#df_fates = pd.DataFrame(fates.X, 
#    columns=fates.var_names if hasattr(fates, 'var_names') else fates.obs.columns,
#    index=fates.obs_names if hasattr(fates, 'obs_names') else fates.index)
#df_fates.to_csv('leishmania_trajectory/results/wot/fate_matrices/fate_matrix.csv', index=True, index_label="cell_id")
