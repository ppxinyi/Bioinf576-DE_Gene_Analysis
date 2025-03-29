import unittest
from unittest.mock import patch, MagicMock
import pandas as pd

class TestRunDESeq2(unittest.TestCase):
    
    @patch('rpy2.robjects.r')
    def test_run_deseq2(self, mock_r):
        # Mocking the R environment and DESeq2 functions
        
        # Mock the DESeq2 result data
        mock_results = MagicMock()
        mock_results.return_value = pd.DataFrame({
            'log2FoldChange': [1.2, -0.5, 0.3],
            'pvalue': [0.01, 0.05, 0.2],
            'padj': [0.02, 0.07, 0.25]
        })
        
        # Mock the DESeq2 function chain
        mock_r.DESeqDataSetFromMatrix = MagicMock()
        mock_r.DESeq = MagicMock()
        mock_r.results = mock_results

        # Sample input data (counts and metadata)
        counts_data = {
            'sample1': [100, 150, 200],
            'sample2': [120, 160, 180],
            'sample3': [110, 140, 190]
        }
        counts_df = pd.DataFrame(counts_data, index=['gene1', 'gene2', 'gene3'])
        
        col_data = pd.DataFrame({
            'condition': ['treated', 'control', 'treated']
        }, index=['sample1', 'sample2', 'sample3'])
        
        # Call the function
        result_df = run_deseq2(counts_df, col_data)

        # Check if the result is a DataFrame
        self.assertIsInstance(result_df, pd.DataFrame)

        # Check if the expected columns are in the result
        self.assertTrue('log2FoldChange' in result_df.columns)
        self.assertTrue('pvalue' in result_df.columns)
        self.assertTrue('padj' in result_df.columns)
        
        # Check the content of the result (example)
        self.assertEqual(result_df['log2FoldChange'][0], 1.2)
        self.assertEqual(result_df['pvalue'][1], 0.05)

class TestCompareConditions(unittest.TestCase):

    def test_compare_conditions(self):
        mock_df = pd.DataFrame({
            'condition': ['A', 'A', 'B', 'B', 'C', 'C'],
            'log2FoldChange': [1.5, -1.2, 0.8, -0.6, 2.1, -2.3],
            'padj': [0.01, 0.04, 0.2, 0.1, 0.03, 0.06]
        })

        # Call the function for two conditions
        result_two = compare_conditions(mock_df, 'condition', ['A', 'B'], pval_threshold=0.05, logfc_threshold=1.5)

        # Call the function for three conditions
        result_three = compare_conditions(mock_df, 'condition', ['A', 'B', 'C'], pval_threshold=0.05, logfc_threshold=1.5)

        # Check if the result is a DataFrame
        self.assertIsInstance(result_two, pd.DataFrame)
        self.assertIsInstance(result_three, pd.DataFrame)

        # Check that the correct number of rows are returned
        self.assertEqual(len(result_two), 1)  
        self.assertEqual(len(result_three), 2)  


if __name__ == '__main__':
    unittest.main()
