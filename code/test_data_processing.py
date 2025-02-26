import unittest
import numpy as np
import pandas as pd
from data_provessing import add

class TestNormalization(unittest.TestCase):

    def setUp(self):
        self.gene_lengths = pd.Series([1000, 2000, 1500], index=["GeneA", "GeneB", "GeneC"])
        self.counts_df = pd.DataFrame(
            {
                "Sample1": [100, 200, 300],
                "Sample2": [400, 500, 600]
            },
            index=["GeneA", "GeneB", "GeneC"]
        )

    def test_raw_normalization(self):
        normalized_df = normalize_counts(self.counts_df, method="raw")
        self.assertFalse(normalized_df.isnull().values.any(), "No NaN values")
        self.assertEqual(normalized_df.shape, self.counts_df.shape, "keep same shape")

    def test_fpkm_normalization(self):
        normalized_df = normalize_counts(self.counts_df, method="FPKM")
        self.assertFalse(normalized_df.isnull().values.any(), "No NaN values")
        self.assertEqual(normalized_df.shape, self.counts_df.shape, "keep same shape")

    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            normalize_counts(self.counts_df, method="INVALID")

if __name__ == "__main__":
    unittest.main()
