from openpyxl import load_workbook
from collections import Counter

file_path = '/home/foundation/data/foundation/data/data/Mitochondrial Proteins-240724.xlsx'

wb = load_workbook(filename=file_path, data_only=True)
sheet = wb['Selected+glucolysis']
column_letter1 = 'C'
column_letter2 = 'D'

# 获取两列数据
gene_symbol_list = [cell.value for cell in sheet[column_letter1][1:]]
gene_list = [cell.value for cell in sheet[column_letter2][1:]]

wb.close()

print(len(gene_list))
print(len(gene_symbol_list))

# # # 找出两列数据的差异
# diff1 = set(gene_symbol_list) - set(gene_list)
# diff2 = set(gene_list) - set(gene_symbol_list)

# print("在 gene_symbol_list 中但不在 gene_list 中的基因：", diff1)
# print("在 gene_list 中但不在 gene_symbol_list 中的基因：", diff2)

# print(len(diff1))
# print(len(diff2))

# # 检查一一对应的关系
# differences = []
# for index, (gene_symbol, gene) in enumerate(zip(gene_symbol_list, gene_list)):
#     if gene_symbol != gene:
#         differences.append((index, gene_symbol, gene))

# print(len(differences))
        
# if differences:
#     print("两列数据不对应的项：")
#     for diff in differences:
#         print(f"第 {diff[0] + 1} 行：gene_symbol_list 为 {diff[1]}，gene_list 为 {diff[2]}")
# else:
#     print("两列数据完全对应")

    
# # 查找重复值及其出现次数
# gene_symbol_counter = Counter(gene_symbol_list)
# gene_list_counter = Counter(gene_list)

# duplicate_genes = {gene: gene_symbol_counter[gene] + gene_list_counter[gene]
#                    for gene in set(gene_symbol_list + gene_list)
#                    if gene_symbol_counter[gene] + gene_list_counter[gene] > 1}

# print("重复的基因及其出现次数：", duplicate_genes)
    
# 合并两列数据，并去除重复
combined_gene_list = list(set(gene_symbol_list + gene_list))

print("合并后的基因列表长度：", len(combined_gene_list))

# # 验证 diff1 和 diff2 中的基因是否都在 combined_gene_list 中
# missing_in_combined = [gene for gene in (diff1 | diff2) if gene not in combined_gene_list]
# if missing_in_combined:
#     print("以下基因在 combined_gene_list 中未找到：", missing_in_combined)
# else:
#     print("所有 diff1 和 diff2 中的基因都在 combined_gene_list 中找到")
    
gene_dict = {gene: index for index, gene in enumerate(combined_gene_list, start=0)}
import json

# 写入 JSON 文件
with open('/home/foundation/data/foundation/data/data/vocab-1162.json', 'w') as json_file:
    json.dump(gene_dict, json_file, indent=4)

print("基因列表已写入 gene_list.json")