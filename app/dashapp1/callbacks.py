import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_html_components as html
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import pandas as pd
import metabolyze as met
from dash.dependencies import Input, Output
from functools import reduce
import dash_html_components as html
import dash_core_components as dcc
import pandas as pd
import metabolyze as met
from functools import reduce
import warnings
import pandas as pd
import itertools
import scipy
import scipy.stats
import numpy as np
from functools import reduce
import re
import numpy 
import subprocess as sp
import os
import sys
import time
import dash
import dash_table
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, Input, State
import plotly.graph_objs as go
from sklearn.preprocessing import scale
import scipy.cluster.hierarchy as shc
import pandas as pd
import os
import matplotlib
matplotlib.use('Agg')  
import itertools
import string
#import RNAseq
import numpy as np
import pandas as pd
from scipy import stats
import numpy as np
import sys, json
from time import sleep
from scipy import *
from scipy import stats
from sklearn.decomposition import PCA
import chart_studio.plotly as py
#%matplotlib inline
import plotly
#plotly.offline.init_notebook_mode() # To embed plots in the output cell of the notebook
import plotly.graph_objs as go
import os
import operator
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from math import sqrt
from itertools import combinations
from matplotlib import rcParams


result = met.Analysis('skeleton_output.tsv','Groups.csv')

def register_callbacks(dashapp):
    
    COLORS20 = [
    '#1f77b4',
    '#aec7e8',
    '#ff7f0e',
    '#ffbb78',
    '#2ca02c',
    '#98df8a',
    '#d62728',
    '#ff9896',
    '#9467bd',
    '#c5b0d5',
    '#8c564b',
    '#c49c94',
    '#e377c2',
    '#f7b6d2',
    '#7f7f7f',
    '#c7c7c7',
    '#bcbd22',
    '#dbdb8d',
    '#17becf',
    '#9edae5',
    ]
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value', 'children'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])
    def update_output(n_clicks,value):
        standard = pd.read_table(result.data)
        detection_column_index = standard.columns.get_loc("detections")

        standard = standard.iloc[:,0:detection_column_index]
        standard.index = standard['Metabolite']
        del standard['Metabolite']


        matrices = []    
        sample_groups = result.get_groups()
        #print (comparison[0])
        #value = value.replace("'", "")

        value_list = value.split(',')
        #value_list_strip = [value.replace(" ", "") for value in value_list]
        value_list_final = [val[1:-1] for val in value_list]
        #value_list = value_list.replace("'", "")
        #test_condition = 'GPI1'
        #ids = (sample_groups[test_condition]) 
        #print('wtf')
        comparison_ids = []

        matrices = []
        for condition in value_list_final:
            #print("condition",condition)
            if condition in sample_groups:
                test = condition
                #real_condition = condition.replace("'", "")
                ids = (sample_groups[test]) 
                #print (ids)
                matrices.append((result.get_imputed_full_matrix(result.get_matrix(ids=ids),param='detected')))
                comparison_ids.append(ids)
        #print(matrices)
        #print("yowtf",matrices[0].shape[1])
        group_sample_number =  int((matrices[0].shape)[1])
        group_sample_number_2 = int(group_sample_number+ ((matrices[1].shape)[1]))

        df_m = reduce(lambda left,right: pd.merge(left,right,left_index=True, right_index=True), matrices)
        intensity_matrix = reduce(lambda left,right: pd.merge(left,right,left_index=True, right_index=True), matrices)

        blank_matrix = pd.DataFrame(result.get_matrix(result.get_ids('Blank')))
        #blank_matrix.to_csv(results_folder+'Tables/'+'blank_intensity.csv')
        blank_threshold = pd.DataFrame(blank_matrix.mean(axis=1)*3)+10000
        blank_threshold['Metabolite'] = blank_threshold.index
        blank_threshold.columns = ['blank_threshold','Metabolite']

        #print(df_m.head())
        df_m['ttest_pval'] = ((scipy.stats.ttest_ind(df_m.iloc[:, :group_sample_number], df_m.iloc[:, group_sample_number:group_sample_number_2], axis=1))[1])
        df_m['1/pvalue'] = float(1)/df_m['ttest_pval']      
        group_1_df = (pd.DataFrame(df_m.iloc[:, :group_sample_number]))
        group_2_df = (pd.DataFrame(df_m.iloc[:, group_sample_number:group_sample_number_2]))
        
        
        
        df_m[value_list_final[0]+'_Mean'] = (group_1_df.mean(axis=1))
        df_m[value_list_final[1]+'_Mean'] = (group_2_df.mean(axis=1))
        
        df_m['Log2FoldChange'] =  np.log2(((group_1_df.mean(axis=1)))/((group_2_df.mean(axis=1))))
        df_m['LogFoldChange'] =  (((group_1_df.mean(axis=1)))/((group_2_df.mean(axis=1))))
        # df_m['Metabolite'] = df_m.index


        final_df_m = pd.merge(standard, df_m, left_index=True, right_index=True)
        final_df_m = pd.merge(final_df_m,blank_threshold,left_index=True, right_index=True)
        # # Add detection column

        for col in blank_matrix.columns:

            final_df_m[col] = blank_matrix[col].values


        final_df_m['combined_mean'] = (final_df_m[value_list_final[0]+'_Mean']+final_df_m[value_list_final[1]+'_Mean'])/2
        final_df_m['impact_score'] = (((2**abs(final_df_m['Log2FoldChange']))*final_df_m['combined_mean'])/final_df_m['ttest_pval'])/1000000
        final_df_m.impact_score = final_df_m.impact_score.round()
        final_df_m['impact_score'] = final_df_m['impact_score'].fillna(0)



        sig_genes = final_df_m.loc[final_df_m['ttest_pval'] < 0.05].index
        

        # data = cluster_plot(intensity_matrix,sig_genes)


        detection_dict = {}
                
        comparison_matrix = group_1_df.join(group_2_df, how='outer')
        
        
        for index, row in comparison_matrix.iterrows():
            test_list = []
            #print (row)
            #print(index)
            row_intensity = (pd.DataFrame(row))
            blankthresh = blank_threshold.loc[index, ['blank_threshold']][0]
            detected = (row_intensity[row_intensity > float(blankthresh)].count())
            detected = (detected[0])
            detection_dict[index] = detected

        detection_df = pd.DataFrame(list(detection_dict.items()))
        detection_df.columns = ['Metabolite','Detection']
        detection_df.index = detection_df['Metabolite']

        final_df_m = pd.merge(final_df_m,detection_df,left_index=True, right_index=True)
        # volcano df
        
        #final_df_m = final_df_m.fillna('NA')
        final_df_m = final_df_m[np.isfinite(final_df_m['ttest_pval'])]

        return (final_df_m.to_json(),intensity_matrix.to_json())


    @dashapp.callback(
        dash.dependencies.Output('output-data-upload', 'children'),
        [dash.dependencies.Input('intermediate-value', 'children')])

    def update_table(jsonified_cleaned_data):
        dff = pd.read_json(jsonified_cleaned_data[0]) # or, more generally json.loads(jsonified_cleaned_data)

        table = html.Div([dash_table.DataTable(
        id = 'output-data-upload',
        columns=[{"name": i, "id": i} for i in dff.columns],
        data=dff.to_dict('records'),
        editable=False,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        #fixed_rows={ 'headers': True, 'data': 0 },
        row_selectable="multi",
        row_deletable=True,
        selected_columns=[],
        selected_rows=[],
        page_action="native",
        style_table={'overflowX': 'scroll'},
        export_format='xlsx',
        export_headers='display',
        merge_duplicate_headers=True)],
        className='twelve columns')
        print(dir(table))
        return(table)


    # @dashapp.callback(
    # Output('volcano_graph', 'children'),
    # [Input('output-data-upload', 'derived_virtual_row_ids'),
    #  Input('output-data-upload', 'selected_row_ids'),
    #  Input('output-data-upload', 'active_cell')])
    
    @dashapp.callback(
        dash.dependencies.Output('volcano_graph', 'children'),
        [dash.dependencies.Input('intermediate-value', 'children')])
    #row_ids, selected_row_ids, active_cell

    # @dashapp.callback(
    #     Output('volcano_graph', 'children'),
    #     [Input('output-data-upload', 'children'),
    #     Input('output-data-upload', 'children'),
    #     Input('output-data-upload', 'derived_virtual_selected_rows')])

    def update_volcano(jsonified_cleaned_data):
        

        final_df_m = pd.read_json(jsonified_cleaned_data[0])
        volcano_df = final_df_m[np.isfinite(final_df_m['ttest_pval'])]
        sig_volcano = volcano_df.loc[volcano_df['ttest_pval'] < 0.05]
        insig_volcano = volcano_df.loc[volcano_df['ttest_pval'] > 0.05]


        volcano_fig_1 = go.Scatter(
        #x = -np.log10(comparison_insig['ttest_pval']),


        x = sig_volcano['Log2FoldChange'],
        y = -np.log10(sig_volcano['ttest_pval']),
        #y = comparison_insig['Log2FoldChange'],
        text = sig_volcano.index,
        name = 'Significant',
        mode = 'markers',
        marker = dict(
            size = 10,
            color = 'rgba(152, 0, 0, .8)',
            line = dict(
                width = 2,)))

        volcano_fig_2 = go.Scatter(
        #x = -np.log10(comparison_insig['ttest_pval']),


        x = insig_volcano['Log2FoldChange'],
        y = -np.log10(insig_volcano['ttest_pval']),
        #y = comparison_insig['Log2FoldChange'],
        text = insig_volcano.index,
        name = 'Insignificant',
        mode = 'markers',
        marker = dict(
            size = 10,
            color = 'rgba(255, 182, 193, .9)',
            line = dict(
                width = 2,)))


        volcano_data = [volcano_fig_1,volcano_fig_2]
        
        volcano = html.Div([
                    #html.H3('Volcano Plot',style={'textAlign': 'center'}),
                    dcc.Graph(id='g2', figure={
                    'data': volcano_data,
                    'layout': {
                        'clickmode': 'event+select','width': '450', 'display': 'inline-block', 'height': '450', 'padding': '0',

                    },
                })], className="four columns")

        return(volcano)


    @dashapp.callback(
        dash.dependencies.Output('heatmap_graph', 'children'),
        [dash.dependencies.Input('intermediate-value', 'children')])

    def update_heatmap(jsonified_cleaned_data):
        final_df_m = pd.read_json(jsonified_cleaned_data[0])
        matrix = pd.read_json(jsonified_cleaned_data[1])
        sig_genes = final_df_m.loc[final_df_m['ttest_pval'] < 0.05].index
        def cluster_plot(matrix,genes):

            values = matrix

            values = values.loc[values.index.isin(genes)]
            values_t = values.T
            values_t.index.name = 'Sample'

            #values_t.to_csv('values_t.csv')
            #plt.figure(figsize=(10, 7))  
            #plt.title("Customer Dendograms")  

            dend_metabolite = shc.dendrogram(shc.linkage(values, method='ward'),labels=values.index)  
            dend_metabolite_order = dend_metabolite['ivl']


            dend_sample = shc.dendrogram(shc.linkage(values_t, method='ward'),labels=values_t.index)  
            dend_sample_order = dend_sample['ivl']

            df = values[dend_sample_order]
            df = df.reindex(dend_metabolite_order)



            df_zscore = df.apply(
                        lambda V: scale(V,axis=0,with_mean=True, with_std=True,copy=False),axis=1)


            trace = go.Heatmap(z=np.array(df_zscore),x=dend_sample_order,y=dend_metabolite_order)
            data=[trace]
            return(data)
        data = cluster_plot(matrix,sig_genes)

        heatmap_figure = html.Div([
            #html.H3('Heatmap',style={'textAlign': 'center'}),
            dcc.Graph(id='g3', figure={
            'data': data,
            'layout': {
                'clickmode': 'event+select','width': '450', 'display': 'inline-block', 'height': '450', 'padding': '0',

            },
        })], className="four columns")

        return(heatmap_figure)

    @dashapp.callback(
        dash.dependencies.Output('pca_graph', 'children'),
        [dash.dependencies.Input('intermediate-value', 'children')],
        [dash.dependencies.State('my-dropdown', 'value')])

    def update_pca(jsonified_cleaned_data,value):
        
        matrix = pd.read_json(jsonified_cleaned_data[1])
       
        value_list = value.split(',')
        #value_list_strip = [value.replace(" ", "") for value in value_list]
        value_list_final = [val[1:-1] for val in value_list]
        def perform_PCA(fpkmMatrix, standardize=3, log=True):
        ## preprocessing of the fpkmMatrix
            if log:
                fpkmMatrix = np.log10(fpkmMatrix + 1.)
            if standardize == 2: # standardize along rows/genes
                fpkmMatrix = stats.zscore(fpkmMatrix, axis=1)
            elif standardize == 1: # standardize along cols/samples
                fpkmMatrix = stats.zscore(fpkmMatrix, axis=0)

            ## remove genes with NaNs
            fpkmMatrix = fpkmMatrix[~np.isnan(np.sum(fpkmMatrix, axis=1))]

            pca = PCA(n_components=None)
            ## get variance captured
            pca.fit(fpkmMatrix.T)
            variance_explained = pca.explained_variance_ratio_[0:3]
            variance_explained *= 100
            ## compute PCA and plot
            pca_transformed = pca.transform(fpkmMatrix.T)
            return (variance_explained, pca_transformed)

        def pca_plot(matrix,samplesheet,conditions):

            # directory = matrix_path.split('/')[0]
            # print(directory)
            expr_df = matrix.astype(float)

            variance_explained, pca_transformed = perform_PCA(expr_df.values, standardize=2, log=True)
            

            meta_df = pd.read_csv(samplesheet)
            meta_df = meta_df.loc[meta_df['Group'] != 'Blank']
            meta_df = meta_df.loc[meta_df['Group'].isin(conditions)]

            
            meta_df['File'].str.replace('.mzXML','')

            meta_df.index = meta_df['id']
            #print(meta_df)
            meta_df['x'] = pca_transformed[:,0]
            meta_df['y'] = pca_transformed[:,1]
            meta_df['z'] = pca_transformed[:,2]

            conditions = meta_df['Group'].unique().tolist()
            #platforms = meta_df['platform'].unique().tolist()
            SYMBOLS = ['circle']
            COLORS = COLORS20

            data = [] 
            #print(meta_df)

            for (condition), meta_df_sub in meta_df.groupby(['Group']):
                #print(condition)
                # Iteratate through samples grouped by condition and platform
                display_name = '%s' % (condition)
                # Initiate a Scatter3d instance for each group of samples specifying their coordinates
                # and displaying attributes including color, shape, size and etc.
                trace = go.Scatter3d(
                    x=meta_df_sub['x'],
                    y=meta_df_sub['y'],
                    z=meta_df_sub['z'],
                    text=meta_df_sub.index,
                    mode='markers',
                    marker=dict(
                        size=10,
                        color=COLORS[conditions.index(condition)], # Color by infection status
                        opacity=.8,
                    ),
                    name=display_name,
                )
                
                data.append(trace)

            
            return(data)
        pca_fig = pca_plot(matrix,result.samplesheet,value_list_final)
        pca_div = html.Div([
                #html.H3('PCA',style={'textAlign': 'center'}),
                dcc.Graph(id='g4', figure={
                'data': pca_fig,
                'layout': {
                    'clickmode': 'event+select','width': '500', 'display': 'inline-block', 'height': '450', 'padding': '0',

                },
            })], className="four columns")

        return(pca_div)