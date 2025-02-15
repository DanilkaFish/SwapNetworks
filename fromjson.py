import json
import numpy as np
dic = {}
help_dict = {"swap_sh": ["green", 'pentagon*'],
            "swap_sh_inv" : ["green", 'pentagon*'],
            "jw" : ["blue", 'triangle*', "JW"],
            "jw_lex" :  ["red", 'triangle*'],
            "bk" :  ["blue", 'diamond*', "BK"],
            "bk_lex" : ["red", 'diamond*'],
            "jw_opt" : ["blue", 'square*', "JW OPT"],
            "jw_opt_lexic" : ["red", 'square*'],
            # "swap_yo" : ["black", 'pentagon*'],
            # "swap_yo_inv" : ["black", 'pentagon*'],
            "swap 2xn" : ["purple", 'pentagon*', "SWAP 2xn"],
            "swap gen yor" : ["green", 'pentagon*', "GSN Yordan"],
            "swap gen short" : ["black", 'pentagon*', "GSN short"],}
with open('data/SimParZ8.json', 'r') as rf:
    data = json.load(rf)
    for obj in data:
        # dic[obj["noise gates"]]
        if obj["name"] in dic:
            dic[obj["name"]].append((obj["energy"], obj["noise gates"]))
        else:
            dic[obj["name"]] = [(obj["energy"], obj["noise gates"])]
            print(obj["name"])
    ref_en = data[0]["ref_ener"]

# print(dic)

def gen_latex(dic, ref_en, probs):
    txt = ""
    for name in dic:
        txt += f"\\addplot[color={help_dict[name][0]},mark={help_dict[name][1]}]\n coordinates\u007b"
        for index, en in enumerate(dic[name]):
            txt += f"({probs[index]:.5}, {abs(en[0]-ref_en)/abs(ref_en):.5}) "
        txt += "};\n"
        txt += f"\\addlegendentry \u007b{help_dict[name][2]} {({dic[name][-1][1]})}\u007d \n"
    print(txt)
probs = np.flip(np.geomspace(0.0002, (0.01), 10))

gen_latex(dic, ref_en, probs)
