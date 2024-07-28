import os

#TODO
def extract_lines(file_path, match_str1, match_str2):
    lines = []
    match_count = 0
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            if match_count == 0 and match_str1 in line:
                lines.append(line.strip())
                match_count += 1
            elif match_count == 1 and match_str2 in line:
                lines.append(line.strip())
                match_count += 1
            if match_count == 2:
                break
    return lines


def process_files(input_dir, output_file, match_str1, match_str2):
    with open(output_file, "w", encoding="utf-8") as out_file:
        for file_name in os.listdir(input_dir):
            if file_name.startswith("888") and file_name.endswith("-1.1.txt"):
                file_path = os.path.join(input_dir, file_name)
                lines = extract_lines(file_path, match_str1, match_str2)
                if len(lines) == 2:
                    out_file.write(f"File Name: {file_name}\n")
                    out_file.write(f"{lines[0]}\n")
                    out_file.write(f"{lines[1]}\n")
                    out_file.write("\n")  # Add empty line to separate files


input_dir = "/home/foundation/program/foundation-new/record/txt/"
output_file = "summary-888-1.1.txt"
match_str1 = "Combination:"
match_str2 = "ARI:"

process_files(input_dir, output_file, match_str1, match_str2)