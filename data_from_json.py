import json
from texttable import Texttable
import latextable

f1 = open('./data/numbers_ideal_6.json')
f2 = open('./data/numbers_noise_6.json')
data_ideal = json.load(f1)["data"]
data_noise = json.load(f2)["data"]
data_ideal = sorted(data_ideal, key=lambda x: x["opt"])
# data_noise = sorted(data_noise, key=lambda x: x["opt"])

def obj_to_row(obj):
    def roun(d:float):

        return "{:.5f}".format(d) + '(+)'
    def orb_ro_str(o):
        if sum(o) == 1:
            return " 4q"
        elif sum(o) == 6:
            return " 8q"
        else:
            return " 8q" 
    row = []
    # print(obj)
    print(roun(obj["hf"][0]))
    # print(obj["mol"][1])
    name = '$H_2$' if sum(obj["mol"][1]) == 2 else "$LiH$"
    row.append(name + orb_ro_str(obj["act_orb"]))
    row.append(obj["opt"])
    # row.append(roun(obj["hf"][0]))
    # row.append(roun(obj["ref_ener"]))
    for i in range(7):
        # row.append(roun(obj["energy"][i]) + f" ({obj["iter_num"][i]})")
        row.append(f"{obj["iter_num"][i]}")
    return row
    
def obj_to_row_noise(obj):
    def roun(d:float):

        return "{:.5f}".format(d) + '(+)'
    def orb_ro_str(o):
        if sum(o) == 1:
            return " 4q"
        elif sum(o) == 6:
            return " 8q"
        else:
            return " 8q" 
    row = []
    # print(obj)
    print(roun(obj["hf"][0]))
    # print(obj["mol"][1])
    name = '$H_2$' if sum(obj["mol"][1]) == 2 else "$LiH$"
    row.append(name + orb_ro_str(obj["act_orb"]))
    row.append(obj["opt"])
    # row.append(roun(obj["hf"][0]))
    # row.append(roun(obj["ref_ener"]))
    for i in range(7):
        # row.append(roun(obj["energy"][i]) + f" ({obj["iter_num"][i]})")
        row.append(roun(obj["energy"][i] - obj["ref_ener"]))
    return row

table_1 = Texttable()
table_1.set_cols_align(["c", "c", "c", "c", "c",  "c", "c", "c", "c"])
# table_1.set_cols_valign(["t", "m", "b"])
table_1.add_row(["Molecule", "Opt","dynamic","jw", "jw_lexic", "bk", "bk_lexic", "jw_opt", "jw_opt_lexic"])
for obj in data_ideal:
    row = obj_to_row(obj)
    table_1.add_row(obj_to_row(obj))

table_2 = Texttable()
table_2.set_cols_align(["c", "c", "c", "c", "c",  "c", "c", "c", "c"])
# table_1.set_cols_valign(["t", "m", "b"])
table_2.add_row(["Molecule", "Opt","dynamic","jw", "jw_lexic", "bk", "bk_lexic", "jw_opt", "jw_opt_lexic"])
for obj in data_noise:
    row = obj_to_row_noise(obj)
    table_2.add_row(row)

print('-- Example 1: Basic --')
print('Texttable Output:')
print(table_1.draw())
print('\nLatextable Output:')
print(latextable.draw_latex(table_1, caption="An example table.", label="table:example_table"))

print('-- Example 2: Basic --')
print('Texttable Output:')
print(table_2.draw())
print('\nLatextable Output:')
print(latextable.draw_latex(table_2, caption="An example table.", label="table:example_table"))
