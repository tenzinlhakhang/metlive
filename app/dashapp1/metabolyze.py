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


warnings.filterwarnings("ignore")

#import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument('-rn', "--rowname", nargs='?', help="Rownames for heatmaps // True or False", const=1, type=str, default='True')
# args = parser.parse_args()

    
class Analysis:
    
    def __init__(self, data,samplesheet):
        
        self.data = 'inputs/'+data
        self.samplesheet = 'inputs/'+samplesheet
#         self.heatmap_rowname = args.rowname


    def input_check(self):
        id_dict = self.get_ids('ID')
        print("Number of Samples:",len(id_dict))
        
        for x,y in id_dict.items():
            print (x,':',y)
        sample_id = self.get_ids('All')
        if len(sample_id) != len(set(sample_id)):
            raise Exception('Error: Check unique Sample IDs in: Groups.csv for error')
        
        skeleton_input = pd.read_table(self.data)
        metabolite_list = skeleton_input['Metabolite']
        if len(metabolite_list) != len(set(metabolite_list)):
            raise Exception('Error: Check Metabolite column for duplicates in : Skeleton_input.tsv')
        
        if self.get_matrix(self.get_ids('All')).isnull().values.any():
            raise Exception('Error: Check for Missing Values in Sample intensities: Skeleton_input.csv')
        
        if len(sample_id) != len(test.get_matrix(test.get_ids('All')).columns):
            raise Exception('Error: Check if Number of Samples in Groups.csv matches Skeleton_input.tsv')
    
        skeleton = self.get_ids('All')
        groups = pd.read_csv(self.samplesheet)['File'].tolist()
        
        if set(groups).issubset(skeleton) == False:
            raise Exception('Samplesheet Sample Names Incorrectly Match Skeleton File Names')
    

        
    def dir_create(self):
        groups = pd.read_csv(self.samplesheet)
        results_folder  = 'DME-results-'+str(len(self.get_ids('True'))) + '-Samples/'
        sub_directories = [results_folder+ subdir for subdir in ['Volcano','Heatmap','Tables','PCA','Inputs','Pathway','Correlation']]
        sub_directories.append(results_folder)
        
        for direc in sub_directories:
            if not os.path.exists(direc):
                os.makedirs(direc)

    
    def get_groups(self):
    # Get corresponding IDs for each group in Groups.csv

        project = pd.read_csv(self.samplesheet)
        grouped_samples = {}

        for condition in (project.Group.unique()):
            if condition != 'Blank':
                test = [x.split('.')[0] for x in project.loc[project['Group'] == condition, 'File'].tolist()]
                grouped_samples[condition] = test
        return (grouped_samples)

    def get_ids(self,full):
        
        # Return sample IDS for all samples including blanks
        if full == 'All':
            skeleton = pd.read_table(self.data)
            
            spike_cols = [col for col in skeleton.columns if 'S' in col]
            spike_cols.pop(0)
            return (list(spike_cols))
        
        # Get all sequence IDS (xml ids) from Groups.csv
        if full == 'True':
            project = pd.read_csv(self.samplesheet)
            project = project.loc[project['Group'] != 'Blank']
            all_samples = [x.split('.')[0] for x in project['File'].tolist()]
            return(all_samples)
        
        if full == 'Sample':
            project = pd.read_csv(self.samplesheet)
            project = project.loc[project['Group'] != 'Blank']
            all_samples = [x.split('.')[0] for x in project['id'].tolist()]
            return(all_samples)
        
        # Get all blank IDS from skeleton output matrix
        if full == 'Blank':
            project = pd.read_csv(self.samplesheet)
            project = project.loc[project['Group'] == 'Blank']
            all_samples = [x.split('.')[0] for x in project['File'].tolist()]
            return (list(all_samples))
        if full == 'ID':
            project = pd.read_csv(self.samplesheet)
            grouped_samples = {}
            
            for condition in (project.id.unique()):

                test = [x.split('.')[0] for x in project.loc[project['id'] == condition, 'File'].tolist()]
                test = ''.join(test)
                grouped_samples[test] = condition
            return(grouped_samples)
    
    def sequence2id(self,result):
        
        ids = self.get_ids('ID')
    
        for x,y in ids.items():
            #print(x,y)
            result.rename(columns={x: y}, inplace=True)
            # Returns matrix based on inputted IDS
        return(result)
    
    def get_matrix(self,ids):
        
        skeleton_outbut_hybrid = pd.read_table(self.data)
        skeleton_outbut_hybrid = skeleton_outbut_hybrid.set_index('Metabolite')
        
        matrix = (skeleton_outbut_hybrid[skeleton_outbut_hybrid.columns.intersection(ids)])
        return (matrix)
    
    def get_imputed_full_matrix(self,full_matrix,param):
        
        blank_matrix = pd.DataFrame(self.get_matrix(self.get_ids('Blank')))
        blank_threshold = pd.DataFrame(blank_matrix.mean(axis=1)*3)+10000
        blank_threshold['Metabolite'] = blank_threshold.index
        blank_threshold.columns = ['blank_threshold','Metabolite']


        test_dictionary = {}
        for index, row in full_matrix.iterrows():
            test_list = []
    #print(index)
            for val in row:
                blankthresh = blank_threshold.loc[index, ['blank_threshold']][0]
                if val < blankthresh:
                    if param == 'detected':
                        test_list.append(blankthresh)
                    if param == 'corrected':
                        test_list.append(0)
                else:
                    test_list.append(val)
            test_dictionary[index] = test_list

        df_test = (pd.DataFrame.from_dict(test_dictionary))
        final = df_test.transpose()
        final.columns = list(full_matrix)
        return(final)

     

    def compile_tests(self,results_folder,full_matrix):
        test_compile = {}


        blank_matrix = pd.DataFrame(self.get_matrix(self.get_ids('Blank')))
        blank_threshold = pd.DataFrame(blank_matrix.mean(axis=1)*3)+10000
        blank_threshold['Metabolite'] = blank_threshold.index
        blank_threshold.columns = ['blank_threshold','Metabolite']

            
            
        for file in os.listdir(results_folder):
            if file.endswith('corrected.csv'):
                #path = os.path.abspath(results_folder+file)
                test = pd.read_csv(results_folder+file,keep_default_na=True)
                test = test.fillna('NA')
                test.index = test['Metabolite']
                columns = ['ttest_pval', 'Log2FoldChange','impact_score']
                changed_names = [file +'_'+ x for x in columns]
                changed_names = [x.replace('.corrected.csv','') for x in changed_names]
                
                df1 = pd.DataFrame(test, columns=columns)
                df1.columns  = changed_names
                test_compile[file] = df1
        
        merged_df = pd.concat(test_compile, axis =1)
        merged_df.columns = [col[1] for col in merged_df.columns]
        
        
        test_dictionary = {}
        for index, row in full_matrix.iterrows():
            test_list = []
        #print(index)
            for val in row:
                blankthresh = blank_threshold.loc[index, ['blank_threshold']][0]
                if val < blankthresh:
                    test_list.append(blankthresh)
                else:
                    test_list.append(val)
            test_dictionary[index] = test_list
            
        df_test = (pd.DataFrame.from_dict(test_dictionary))
        final = df_test.transpose()
        final.columns = list(full_matrix)

            
        detection_dict = {}
        for index, row in final.iterrows():
            test_list = []
            #print (row)
            #print(index)
            row_intensity = (pd.DataFrame(row))
            blankthresh = blank_threshold.loc[index, ['blank_threshold']][0]
            detected = (row_intensity[row_intensity > float(blankthresh)].count())
            detected = (detected[0])
            detection_dict[index] = detected
            
        
        test_dictionary = {}
        for index, row in full_matrix.iterrows():
            test_list = []
        #print(index)
            for val in row:
                blankthresh = blank_threshold.loc[index, ['blank_threshold']][0]
                if val < blankthresh:
                    test_list.append('-')
                else:
                    test_list.append(val)
            test_dictionary[index] = test_list
            
        df_test = (pd.DataFrame.from_dict(test_dictionary))
        new_final = df_test.transpose()
        new_final.columns = list(full_matrix)

        detection_df = pd.DataFrame(list(detection_dict.items()))
        detection_df.columns = ['Metabolite','Detection']
        detection_df.index = detection_df['Metabolite']
        
        #detection_df.to_csv()
#       

        compiled = new_final.join(merged_df, how='outer')
        compiled_final = compiled.join(detection_df, how='outer')

        #passing_df = detection_df.drop('Detection', 1)
    
        return(compiled_final,final)

    def dme_comparisons(self):
        
        sample_groups = self.get_groups()
        groups = pd.read_csv(self.samplesheet)
        unique_groups = [x for x in groups.Group.unique() if x != 'Blank']
        unique_comparisons = []
        
        for L in range(0, len(unique_groups)+1):
            for subset in itertools.combinations(unique_groups, L):
                if len(subset)== 2:
                    unique_comparisons.append(subset)
        

        reversed_groups = []
        for comparison in unique_comparisons:
            reversed_comparison = (tuple(((reversed(comparison)))))
            #print(reversed_comparison)
            reversed_groups.append(reversed_comparison)
        #     print(comparison)
        #     print(reversed_comparison)
        #     print("\n")


        unique_comparisons = unique_comparisons + reversed_groups
        
        return(unique_comparisons)

    
    def t_test(self):
        print("\n")
        print("################")
        print("Pipeline executed:")
      
        self.input_check()
        print("\n")
        print("Creating Directories...")
        print("\n")
        # Create all necessary directories
        self.dir_create()
        
        groups = pd.read_csv(self.samplesheet)
        unique_groups = [x for x in groups.Group.unique()]
               
        # get all unique comparisons from Groups.csv
        unique_comparisons = self.dme_comparisons()

        #Meta Data on Metabolites
        standard = pd.read_table(self.data)
        detection_column_index = standard.columns.get_loc("detections")
        standard = standard.iloc[:,0:detection_column_index]

        # Set directory for results folder 
        results_folder  = 'DME-results-'+str(len(self.get_ids('True'))) + '-Samples/'
        
        
        # Get full matrix of intensity values with Sequence IDS replaced with ID from Groups.csv
        full_matrix = self.get_matrix(self.get_ids(full='True'))
        full_matrix = self.sequence2id(full_matrix)
        full_matrix_name = results_folder+'Tables/'+'Intensity.values.csv'
        detected_matrix_name = results_folder+'Tables/'+'Intensity.detected.values.csv'
        full_matrix.to_csv(full_matrix_name)
        
        corrected_matrix = self.sequence2id(self.get_imputed_full_matrix(self.get_matrix(ids=self.get_ids('True')),param        ='corrected'))
        corrected_matrix.index.name = 'Metabolite'
        corrected_matrix.to_csv(results_folder+'Tables/'+'Intensity.corrected.values.csv')
        
        
        for comparison in unique_comparisons:
            matrices = []    
            sample_groups = self.get_groups()
            #print (comparison[0])
            
            comparison_ids = []
            for condition in comparison:   
                if condition in sample_groups:
                    ids = (sample_groups[condition]) 
                    #print (ids)
                    matrices.append((self.get_imputed_full_matrix(self.get_matrix(ids=ids),param='detected')))
                    comparison_ids.append(ids)
            
            
            sample_ids = [item for sublist in comparison_ids for item in sublist]
            #generate samplesheet just for comparison
            
            
            samplesheet = pd.read_csv(self.samplesheet)

            samplesheet_comparison = samplesheet.loc[samplesheet['File'].isin(sample_ids)]
            
            samplesheet_comparison_name = results_folder+'PCA/samplesheet.csv'
            samplesheet_comparison.to_csv(samplesheet_comparison_name)
            
            #print ((matrices.shape())
            group_sample_number =  int((matrices[0].shape)[1])
            group_sample_number_2 = int(group_sample_number+ ((matrices[1].shape)[1]))
            
            #print(comparison_ids)
            
            pca_matrix =  reduce(lambda left,right: pd.merge(left,right,left_index=True, right_index=True), matrices)
            #pca_matrix = pd.DataFrame(pca_matrix).set_index('Metabolite')
            pca_matrix.index.name = 'Metabolite'
            comparison_pca_name = (results_folder+'PCA/'+comparison[0]+'_vs_'+comparison[1]+'_PCA.html').replace(" ", "")
            comparison_pca = results_folder+'PCA/PCA_matrix.csv'
            
            
            pca_matrix.to_csv(comparison_pca)
            
            proc = sp.Popen(['python','-W ignore','pca.py',comparison_pca,samplesheet_comparison_name,comparison_pca_name])
            matrices.append(pd.DataFrame(self.get_matrix(self.get_ids(full='Blank'))))
            df_m = reduce(lambda left,right: pd.merge(left,right,left_index=True, right_index=True), matrices)
#             print(df_m.head())                  
#              df_blankless = df_m.copy()
            
            #print(group_sample_number,group_sample_number_2)
           # print(df_blankless.head())
            
            #return(df_blankless)
            
            ### Calculate Pearson Correlation 


            def get_correlation(matrix,group):

                temp_pearson_dict ={}
                cov = samplesheet.loc[samplesheet['Group'] == group]['Covariate']

                for row in matrix.iterrows():
                    index, data = row

                    pearson_correl = np.corrcoef(data, cov)[0, 1]
                    temp_pearson_dict[index] = pearson_correl

                pearson_df = pd.DataFrame([temp_pearson_dict]).T
                pearson_df.columns = [group]
                return(pearson_df)
            
            blank_matrix = pd.DataFrame(self.get_matrix(self.get_ids('Blank')))
            blank_matrix.to_csv(results_folder+'Tables/'+'blank_intensity.csv')
            blank_threshold = pd.DataFrame(blank_matrix.mean(axis=1)*3)+10000
            blank_threshold['Metabolite'] = blank_threshold.index
            blank_threshold.columns = ['blank_threshold','Metabolite']

            
            df_m['ttest_pval'] = ((scipy.stats.ttest_ind(df_m.iloc[:, :group_sample_number], df_m.iloc[:, group_sample_number:group_sample_number_2], axis=1))[1])
            df_m['1/pvalue'] = float(1)/df_m['ttest_pval']      
            group_1_df = (pd.DataFrame(df_m.iloc[:, :group_sample_number]))
            group_2_df = (pd.DataFrame(df_m.iloc[:, group_sample_number:group_sample_number_2]))
            
            
            
            
            df_m[comparison[0]+'_Mean'] = (group_1_df.mean(axis=1))
            df_m[comparison[1]+'_Mean'] = (group_2_df.mean(axis=1))
            
            df_m['Log2FoldChange'] =  np.log2(((group_1_df.mean(axis=1)))/((group_2_df.mean(axis=1))))
            df_m['LogFoldChange'] =  (((group_1_df.mean(axis=1)))/((group_2_df.mean(axis=1))))
            
            final_df_m = pd.merge(standard, df_m, on='Metabolite')
            final_df_m = pd.merge(final_df_m,blank_threshold,on='Metabolite')
            # Add detection column

            for col in blank_matrix.columns:

                final_df_m[col] = blank_matrix[col].values
          
            comparison_name = (results_folder+'Tables/'+comparison[0]+'_vs_'+comparison[1]+'.corrected.csv').replace(" ", "")
            
            
            
            
            final_df_m = self.sequence2id(final_df_m)
            
            final_df_m['combined_mean'] = (final_df_m[comparison[0]+'_Mean']+final_df_m[comparison[1]+'_Mean'])/2
            final_df_m['impact_score'] = (((2**abs(final_df_m['Log2FoldChange']))*final_df_m['combined_mean'])/final_df_m['ttest_pval'])/1000000
            final_df_m.impact_score = final_df_m.impact_score.round()
            final_df_m['impact_score'] = final_df_m['impact_score'].fillna(0)

            
            
            ####Calculate Detection
            

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

            final_df_m = pd.merge(final_df_m,detection_df,on='Metabolite')
            
            # Add impact score
            
            
            
            print("Analysis",":",comparison[0]+'_vs_'+comparison[1])
            print('Results Generated: %s'%comparison_name)
            final_df_m = final_df_m.fillna('NA')
            
#             final_df_m = pd.merge(final_df_m,merged_pearson,on='Metabolite',how='outer')
            final_df_m.to_csv(comparison_name)  
            
            
            test = pd.read_csv(comparison_name)
            
            print("Significant Metabolites P-value < 0.05:",len(test.loc[test['ttest_pval'] < 0.05]))
        
            #Generate Volcano
            print("Generating Volcano Plot: %s" %comparison_name)
            proc = sp.Popen(['Rscript','scripts/volcano.plot.R',comparison_name])
           
            
            # Generate heatmaps
            pvalues = [str(0.05)]
            print("Generating Pvalue < 0.05 Heatmap: %s"%comparison_name)
            for pvalue in pvalues:    
            
                proc = sp.Popen(['Rscript','scripts/heatmap.R',comparison_name,pvalue,'TRUE'])
             
            # Generate heatmap with all expressed metabolites


            print("\n")
        
            # Generate 3-D PCA
            
        print("Compiling Comparison - Results - output: dme.compiled.csv")
        
        compiled, imputed_intensities = self.compile_tests(results_folder+'Tables/',full_matrix)
        compiled = compiled.fillna('-')
        
        
        def change_column_order(df, col_name, index):
            cols = df.columns.tolist()
            cols.remove(col_name)
            cols.insert(index, col_name)
            return df[cols]
            
        compiled.to_csv(results_folder+'Tables/'+'dme.compiled.csv')
        
        dme_meta_data = standard[['Metabolite','Formula','Polarity (z)','mz','ppm','RT','RT_range']]
        dme_meta_data.index = dme_meta_data['Metabolite']
        compiled = pd.merge(dme_meta_data,compiled,on='Metabolite')
        compiled = change_column_order(compiled, 'Detection', 7)
        
        compiled.to_csv(results_folder+'Tables/'+'dme.compiled.csv')
        
        imputed_intensities.index.name = "Metabolite"
        #imputed_intensities = imputed_intensities.rename(columns={ imputed_intensities.columns[0]: "Metabolite" })
    
        imputed_intensities.to_csv(results_folder+'Tables/'+'Intensity.detected.values.csv')
        print("Generating Full Heatmap")
        proc = sp.Popen(['Rscript','scripts/heatmap.full.R',full_matrix_name,'nonimputed'])
        proc = sp.Popen(['Rscript','scripts/heatmap.full.R',detected_matrix_name,'imputed'])
        proc = sp.Popen(['python','-W ignore','pca.py',detected_matrix_name,self.samplesheet,(results_folder+'PCA/'+'PCA.full.html')])

        os.remove(comparison_pca)
        os.remove(samplesheet_comparison_name)
        
        from shutil import copyfile
        
        copyfile('inputs/Groups.csv', results_folder+'Inputs/'+'Groups.csv')
        copyfile('inputs/skeleton_output.tsv', results_folder+'Inputs/'+'skeleton_output.tsv')

        table_directory = results_folder+'Tables'
        print("resultsfolder path")
        

        print('#######')
#         for file in os.listdir(results_folder+'Tables'):
#             if file.endswith('corrected.csv'):
        path = os.path.abspath(results_folder+'Tables')
        output_path = os.path.abspath(results_folder+'Pathway')

        proc = sp.Popen(['Rscript','scripts/pathway.R',path,output_path])
#                 time.sleep(2)
        impact_folder = results_folder + 'Tables/dme.compiled.csv'
        proc = sp.Popen(['python','scripts/impact.correlation.py', impact_folder])


        proc = sp.Popen(['python','scripts/sig.genes.py',path])
        print("\n")
        print("\n")
        print("\n")
        print("#######")
        print("\n")
        print("\n")
        print("\n")


