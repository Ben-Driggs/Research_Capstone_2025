import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares, lsq_linear


def graph_scatterplot(title, x_label, x_axis, y_label, y_axis, color, count):
    plt.subplot(4, 2, count)
    plt.scatter(x_axis, y_axis, color=color)
    # plt.title(title)
    plt.xlim(0, 5)
    plt.ylim(-1, 10)
    plt.xlabel(x_label, labelpad=5)
    plt.ylabel(y_label, labelpad=5)
    return plt


def main():
    mods = ['m', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
    args = sys.argv[1:]
    results = {}
    num_filtered = {}
    conditions = {}
    
    fig = plt.figure(figsize=(10, 10))
    std_distances = {}
    
    # mode = "compare"
    mode = "calc"
    
    count = 1
    for d in args:
        
        d = d.split(',')
        group = d[0][:6]
        # group = "Diet_A"
        
        # set plot colors
        if group == "Diet_A":
            title = "LP AL"
            color = "grey"
            color2 = "blue"
            # color3 = "navy"
            # color4 = "deepskyblue"
        elif group == "Diet_C":
            title = "LP CR"
            color = "grey"
            color2 = "red"
            # color3 = "darkred"
            # color4 = "darkorchid"
        elif group == "Diet_F":
            title = "HP AL"
            color = "grey"
            color2 = "green"
            # color3 = "darkgreen"
            # color4 = "mediumseagreen"
        else:
            title = "HP CR"
            color = "grey"
            color2 = "orange"
            # color3 = "crimson"
            # color4 = "mediumvioletred"
        
        emp_df = pd.read_csv(d[0], sep='\t', low_memory=False)
        lit_df = pd.read_csv(d[1], sep='\t', low_memory=False)
        
        emp_df = emp_df[emp_df.loc[:, "time"] > 25.0]
        lit_df = lit_df[lit_df.loc[:, "time"] > 25.0]
        
        columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'n_value', 'abundances', 'n_value_stddev']
        l_columns = ['Protein ID', 'Sequence', 'cf', 'mz', 'abundances', 'n_value']
        # columns = ['Protein ID', 'sequence', 'Elemental Composition', 'n_value', 'stddev']
        # l_columns = ['Protein ID', 'sequence', 'Elemental Composition', 'number of possible labeling sites', 'stddev']
        emp_df = emp_df[columns]
        lit_df = lit_df[l_columns]
        
        # merge dataframes
        emp_df['id'] = emp_df['Protein ID'] + emp_df['Sequence'] + str(emp_df['mz']) + emp_df['cf']
        lit_df['id'] = lit_df['Protein ID'] + lit_df['Sequence'] + str(lit_df['mz']) + lit_df['cf']
        emp_df = emp_df.drop(columns=['Protein ID', 'cf', 'mz'])
        lit_df = lit_df.drop(columns=['Protein ID', 'Sequence', 'cf', 'mz', 'abundances'])
        merged_df = pd.merge(emp_df, lit_df, how="outer", on=['id'], suffixes=('_emp', '_lit'))
        
        merged_df = merged_df.drop_duplicates()
        
        merged_df = merged_df[merged_df['n_value_emp'] != "no valid time points"]
        
        # remove any peptides where the n-value is larger than the number of possible labeled sites
        aa_df = pd.read_csv("aa_labeling_sites.tsv", sep='\t')
        max_sites_dict = aa_df.iloc[2, 1:].to_dict()

        def sum_labeling_sites(seq):
            return sum(max_sites_dict.get(a, 0) for a in seq)

        merged_df['max_label_sites'] = merged_df['Sequence'].apply(sum_labeling_sites)
        
        merged_df['n_value_stddev'] = merged_df['n_value_stddev'].astype(float)
        merged_df['n_value_emp'] = merged_df['n_value_emp'].astype(float)
        merged_df['max_label_sites'] = merged_df['max_label_sites'].astype(float)
        merged_df = merged_df[merged_df['n_value_emp'] <= merged_df['max_label_sites']]
        
        # remove standard deviations of zero?
        # merged_df = merged_df[merged_df['n_value_emp'] < (0.9 * merged_df['max_label_sites'])]
        
        # remove any rows with modified peptides
        e_mask = merged_df['Sequence'].str.contains('|'.join(mods))
        merged_df = merged_df[~e_mask]
        
        # filter by abundance
        # sum abundances and take top 50%
        merged_df.loc[:, 'sum_abundances'] = merged_df['abundances'].apply(
            lambda v: sum([float(value) for value in v[1:-1].split(', ')]))

        emp_threshold = merged_df['sum_abundances'].quantile(0.25)

        merged_df = merged_df[merged_df['sum_abundances'] >= emp_threshold]
        
        # emp_threshold = merged_df['sum_abundances'].quantile(0.75)
        # merged_df = merged_df[merged_df['sum_abundances'] <= emp_threshold]
        
        # filter out noise by setting a limit on n_value standard deviation
        merged_df = merged_df[merged_df['n_value_lit'] != "no valid time points"]
        merged_df = merged_df[merged_df['n_value_emp'] != "no valid time points"]
        filtered_df = merged_df[merged_df.loc[:, 'n_value_stddev'] <= 0.05]
        # filtered_df = merged_df.copy()
        
        cols = ['n_value_emp', 'n_value_lit']
        for c in cols:
            merged_df[c] = merged_df[c].astype(float).round(2)
            filtered_df[c] = filtered_df[c].astype(float).round(2)
        
        x = merged_df['n_value_lit']
        y = merged_df['n_value_emp']
        f_x = filtered_df['n_value_lit']
        f_y = filtered_df['n_value_emp']
        
        if mode == "compare":
            plt.subplot(2, 2, count)
            plt.scatter(x, y, color=color, alpha=0.50, label='Unfiltered')
            coeffs1 = np.polyfit(x, y, 1)
            fit1 = np.poly1d(coeffs1)
            plt.plot(x, fit1(x), color="darkgrey", alpha=0.5, linewidth=2)
    
            distances = np.abs(y - x) / np.sqrt(2)
            std_distance = np.std(distances)
            std_distances[group] = std_distance
    
            d_mask = distances <= std_distance
            filt_x = x[d_mask]
            filt_y = y[d_mask]
    
            # plt.scatter(filt_x, filt_y, color=color2, alpha=0.2, label='Filtered')
            # coeffs2 = np.polyfit(filt_x, filt_y, 1)
            # fit2 = np.poly1d(coeffs2)
            # plt.plot(filt_x, fit2(filt_x), color=color2, linewidth=2, label="Filtered Linear Fit")
    
            plt.scatter(f_x, f_y, color=color2, alpha=0.5, label='Filtered')
            coeffs2 = np.polyfit(f_x, f_y, 1)
            fit2 = np.poly1d(coeffs2)
            plt.plot(f_x, fit2(f_x), color=color2, linewidth=2)
    
            upper = [xs + (std_distance * np.sqrt(2)) for xs in list(range(0, 100))]
            lower = [xs - (std_distance * np.sqrt(2)) for xs in list(range(0, 100))]
    
            # plt.plot([0, 100], [upper[0], upper[99]], color=color2, linestyle='--', alpha=0.5)
            # plt.plot([0, 100], [lower[0], lower[99]], color=color2, linestyle='--', alpha=0.5)
            plt.plot([0, 100], [0, 100], color='black', alpha=0.75, linestyle='--')
    
            # plt.xlim(0, 100)
            # plt.ylim(0, 100)
            # plt.title(f"{group}")
            # plt.xlabel(f"{title} Peptide Literature N-Values")
            # plt.ylabel(f"{title} Peptide Empirical N-Values")
            plt.legend(title=title)
            plt.xlim(0, 100)
            plt.ylim(0, 250)
            count += 1
            print(len(f_x))
            
        else:
            no_mods = True
            if no_mods:
                amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                               'Y']
            else:
                amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',
                               'Y',
                               'm', 'q', '1', '2', '3', '4', 's', 't', 'y', 'k', 'c', 'o']
            
            # set plot colors
            if group == "Diet_A":
                color = "blue"
            elif group == "Diet_C":
                color = "red"
            elif group == "Diet_F":
                color = "green"
            else:
                color = "orange"
                
            # create numpy arrays with sequences as rows, and amino acid counts as columns
            emp_aa_matrix = np.zeros((len(f_x), len(amino_acids)), dtype=int)
            lit_aa_matrix = np.zeros((len(f_y), len(amino_acids)), dtype=int)
            
            for i, peptide in enumerate(filtered_df['Sequence'].values):
                for aa in peptide:
                    if aa in amino_acids:
                        emp_aa_matrix[i, amino_acids.index(aa)] += 1
            
            for i, peptide in enumerate(filtered_df['Sequence'].values):
                for aa in peptide:
                    if aa in amino_acids:
                        lit_aa_matrix[i, amino_acids.index(aa)] += 1
            
            emp_n_values = np.array(f_y, dtype=float)
            lit_n_values = np.array(f_x, dtype=float)
            
            # https://numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html#numpy-linalg-lstsq
            # emp_amino_acid_values, emp_residuals, emp_rank, emp_s = np.linalg.lstsq(emp_aa_matrix, emp_n_values, rcond=None)
            # print("Empirical Amino Acid Values:", emp_amino_acid_values)
            
            result1 = lsq_linear(lit_aa_matrix, lit_n_values)
            lit_amino_acid_values = result1.x
            # lit_amino_acid_values, lit_residuals, lit_rank, lit_s = np.linalg.lstsq(lit_aa_matrix, lit_n_values, rcond=None)
            print("Literature Amino Acid Values:", lit_amino_acid_values)
            
            aa_nv = list(aa_df.iloc[0, :].values)
            
            x0 = np.array(aa_nv[2:-12])
            
            def residuals(xs, A, b):
                return np.dot(A, xs) - b
            
            # set bounds for least squares solution
            # lower bound is 0, upper bound is determined by number of possile labeling locations
            bounds = [
                np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
                np.array([4, 3, 3, 5, 8, 2, 5, 10, 9, 10, 8, 3, 7, 5, 7, 3, 5, 8, 8, 7])
            ]
            
            # result = least_squares(residuals, x0, args=(emp_aa_matrix, emp_n_values))
            result = lsq_linear(emp_aa_matrix, emp_n_values, bounds=bounds)
            
            # run monte carlo simulation to try and improve results
            num_samples = 1000
            solutions = []
            for _ in range(num_samples):
                sample = emp_n_values + np.random.uniform(-0.25, 0.25, size=emp_n_values.shape)
                b_sample = emp_n_values + sample
                
                # adjust bounds to account for added noise
                bounds[0] -= sample
                bounds[1] += sample
                
                solution = lsq_linear(emp_aa_matrix, b_sample, bounds=bounds)
                solutions.append(solution.x)
                
            solutions = np.array(solutions)
            
            # calculate mean and standard deviation
            mc_amino_acid_values = np.mean(solutions, axis=0)
            mc_stddev = np.std(solutions, axis=0)
            mc_stddev = [round(x, 2) for x in mc_stddev]
            
            print("Monte Carlo means:", mc_amino_acid_values)
            print("Monte Carlo std devs:", mc_stddev)

            emp_amino_acid_values = result.x
            # print('Empirical Amino Acid Values:', emp_amino_acid_values)
            num_filtered[group] = len(emp_n_values)
            conditions[group] = round(np.linalg.cond(emp_aa_matrix), 2)
            
            if no_mods:
                lit_graph = graph_scatterplot(f"{group} Literature N-Values", "Literature N-Values", aa_nv[2:-12],
                                              "Estimated N-Values", lit_amino_acid_values, color, count)
                count += 1
                emp_graph = graph_scatterplot(f"{group} Empirical N-Values", "Literature N-Values", aa_nv[2:-12],
                                              "Empirical N-Values", emp_amino_acid_values, color, count)
                monte_carlo_graph = graph_scatterplot(f"{group} Monte Carlo N-Values",
                                                      "Literature N-Values", aa_nv[2:-12], "Monte Carlo N-Values",
                                                      mc_amino_acid_values, color, count)
            else:
                lit_graph = graph_scatterplot(f"{group} Literature N-Values", "Literature N-Values", aa_nv[2:],
                                              "Estimated N-Values", lit_amino_acid_values, color, count)
                count += 1
                emp_graph = graph_scatterplot(f"{group} Empirical N-Values", "Literature N-Values", aa_nv[2:],
                                              "Empirical N-Values", emp_amino_acid_values, color, count)
                monte_carlo_graph = graph_scatterplot(f"{group} Monte Carlo N-Values, std={mc_stddev}",
                                                      "Literature N-Values", aa_nv[2:-12], "Monte Carlo N-Values",
                                                      mc_amino_acid_values, color, count)
            count += 1
            results[f"{group}_empirical_n_value"] = [round(v, 2) for v in emp_amino_acid_values]
            results[f"{group}_mc_n_values"] = [round(v, 2) for v in mc_amino_acid_values]
            results[f"{group}_mc_stddevs"] = mc_stddev
            if group == "Diet_G":
                results["literature_n_value"] = [round(v, 2) for v in lit_amino_acid_values]
        
    if mode == "compare":
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plt.tight_layout()
        plt.show()
    else:
        plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
        plt.tight_layout(h_pad=12, w_pad=6)
        plt.show()
        amino_acid_names = [
            'Alanine',
            'Cysteine',
            'Aspartic Acid',
            'Glutamic Acid',
            'Phenylalanine',
            'Glycine',
            'Histidine',
            'Isoleucine',
            'Lysine',
            'Leucine',
            'Methionine',
            'Asparagine',
            'Proline',
            'Glutamine',
            'Arginine',
            'Serine',
            'Threonine',
            'Valine',
            'Tryptophan',
            'Tyrosine'
        ]
        
        df = pd.DataFrame(data=results, index=amino_acid_names)
        df.to_csv(path_or_buf="table_data/amino_acid_n_values", sep='\t')
    
    print("Total:", num_filtered)
    print("Conditions:", conditions)
    sys.exit()


if __name__ == "__main__":
    main()

# TODO: could we implement residuals for amino acids that we don't expect to change much?
# TODO: what is the range of noise we should use for Monte Carlo?
# filter out high standard deviations
# make tablular form of the AA n-values
# recalculate peptide n-values from empirical AA n-values
# filter top 50% abundant peptides
