import scanpy as sc
from anndata import AnnData
import pandas as pd
from itertools import product
# change TODO
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
        # 分割列名字符串并去除首尾空格
        column_names = [col.strip() for col in columns.split(',')]
        
        # 遍历所有可能的组合
        unique_combinations = [adata.obs[col].unique() for col in column_names]
        
        for combination in product(*unique_combinations):
            # 构建条件筛选
            condition = pd.Series([True] * adata.shape[0], index=adata.obs.index)
            for col, value in zip(column_names, combination):
                condition &= (adata.obs[col] == value)
            
            # 过滤数据
            filtered_adata = adata[condition, :]
            
            # 检查过滤后的数据大小
            if 5000 <= filtered_adata.shape[0] <= 20000:
                print("A batch was successfully extracted")
                # 打印这个时候组合的值
                print(f"Combination: {combination}")

                    
                return filtered_adata

        # 如果没有找到符合条件的子集
        print("No batch found with cell count between 5000 and 20000")
        return None

    except Exception as e:
        print(f"Failed to filter data: {e}")
        return None