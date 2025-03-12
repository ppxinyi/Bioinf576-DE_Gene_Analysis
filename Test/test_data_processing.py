import unittest
import pandas as pd
import numpy as np
from data_provessing import normalize_counts, compute_z_scores, filter_low_variance_genes


class TestProcessingFunctions(unittest.TestCase):
    def setUp(self):
        # Sample gene expression data
        self.counts_df = pd.DataFrame({
            "Sample1": [100, 200, 300, 0, 500],
            "Sample2": [50, 100, 150, 0, 250],
            "Sample3": [80, 160, 240, 0, 400]
        }, index=["Gene1", "Gene2", "Gene3", "Gene4", "Gene5"])
        
        # Sample gene lengths
        self.gene_lengths = pd.Series({
            "Gene1": 1000,
            "Gene2": 1500,
            "Gene3": 2000,
            "Gene4": 500,
            "Gene5": 1200
        })

    # Test normalize_counts with "raw" method
    def test_normalize_counts_raw(self):
        normalized_df = normalize_counts(self.counts_df, self.gene_lengths, method="raw")
        self.assertIsInstance(normalized_df, pd.DataFrame)
        self.assertEqual(normalized_df.shape, self.counts_df.shape)  # Ensure dimensions remain unchanged
        self.assertFalse(normalized_df.isnull().values.any())  # Ensure no NaN values

    # Test normalize_counts with "FPKM" method
    def test_normalize_counts_fpkm(self):
        normalized_df = normalize_counts(self.counts_df, self.gene_lengths, method="FPKM")
        self.assertIsInstance(normalized_df, pd.DataFrame)
        self.assertEqual(normalized_df.shape, self.counts_df.shape)
        self.assertFalse(normalized_df.isnull().values.any())

    # Test normalize_counts with an invalid method
    def test_normalize_counts_invalid_method(self):
        with self.assertRaises(ValueError):
            normalize_counts(self.counts_df, self.gene_lengths, method="INVALID")

    # Test compute_z_scores
    def test_compute_z_scores(self):
        z_scores = compute_z_scores(self.counts_df)
        self.assertIsInstance(z_scores, pd.DataFrame)
        self.assertEqual(z_scores.shape, self.counts_df.shape)
        # Check that mean of each row is close to 0 and std is close to 1
        self.assertTrue(np.allclose(z_scores.mean(axis=1), 0, atol=1e-6))
        self.assertTrue(np.allclose(z_scores.std(axis=1, ddof=0), 1, atol=1e-6))

    # Test filter_low_variance_genes
    def test_filter_low_variance_genes(self):
        filtered_df = filter_low_variance_genes(self.counts_df, threshold=1000)  # High threshold to remove all genes
        self.assertIsInstance(filtered_df, pd.DataFrame)
        self.assertLess(filtered_df.shape[0], self.counts_df.shape[0])  # Ensure some genes were removed

if __name__ == "__main__":
    unittest.main()
