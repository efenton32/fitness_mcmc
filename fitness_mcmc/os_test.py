import os

file_name = "test.txt"
start = os.path.abspath(__file__)
start_dir = os.path.dirname(start)
print(start_dir)

sub_dir = "raw_data"
data_path = os.path.join(start_dir, sub_dir, file_name)
print(data_path)

f = open(data_path, "w")
f.write("success")
f.close()