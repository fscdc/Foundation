import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product

def extract_batch(adata: AnnData, columns: str) -> AnnData:
    """
    Filter data based on the specified columns' values.

    Parameters:
    adata (AnnData): Input AnnData object.
    columns (str): Comma-separated column names to filter by.

    Returns:
    AnnData: Filtered AnnData object.
    """
    try:
        # Split column names string and strip leading/trailing whitespace
        column_names = [col.strip() for col in columns.split(',')]

        # Iterate over all possible combinations
        unique_combinations = [adata.obs[col].unique() for col in column_names]

        for combination in product(*unique_combinations):
            # Construct filtering condition
            condition = pd.Series([True] * adata.shape[0], index=adata.obs.index)
            for col, value in zip(column_names, combination):
                condition &= (adata.obs[col] == value)

            # Filter data
            filtered_adata = adata[condition, :]

            # Check the size of the filtered data
            if 5000 <= filtered_adata.shape[0] <= 20000:
                print("A batch was successfully extracted")
                # Print the combination values
                print(f"Combination: {combination}")

                return filtered_adata

        # If no subset found with cell count between 5000 and 20000
        print("No batch found with cell count between 5000 and 20000")
        return None

    except Exception as e:
        print(f"Failed to filter data: {e}")
        return None


def extract_batch2(adata: AnnData, columns: str) -> AnnData:
    """
    Filter data based on the specified columns' values.

    Parameters:
    adata (AnnData): Input AnnData object.
    columns (str): Comma-separated column names to filter by.

    Returns:
    AnnData: Filtered AnnData object.
    """
    try:
        # Split column names string and strip leading/trailing whitespace
        column_names = [col.strip() for col in columns.split(',')]
        
        # Get unique combinations for each column
        unique_combinations = []
        for col in column_names:
            unique_values = adata.obs[col].unique()
            # If column is 'Tissue', exclude 'Normal'
            if col == 'Tissue':
                unique_values = [val for val in unique_values if val != 'Normal']
            unique_combinations.append(unique_values)
        
        # Iterate over all possible combinations
        for combination in product(*unique_combinations):
            # Construct filtering condition
            condition = pd.Series([True] * adata.shape[0], index=adata.obs.index)
            for col, value in zip(column_names, combination):
                condition &= (adata.obs[col] == value)
            
            # Filter data
            filtered_adata = adata[condition, :]
            
            # Check the size of the filtered data
            if 5000 <= filtered_adata.shape[0] <= 20000:
                print("A batch was successfully extracted")
                # Print the combination values
                print(f"Combination: {combination}")

                return filtered_adata

        # If no subset found with cell count between 5000 and 20000
        print("No batch found with cell count between 5000 and 20000")
        return None

    except Exception as e:
        print(f"Failed to filter data: {e}")
        return None
