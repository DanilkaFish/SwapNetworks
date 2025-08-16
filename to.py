import stim
from typing import List, Tuple

import stim._stim_sse2



# ---------- строим U и схему ----------
def clifford_from_sets(in_set: List[str], out_set: List[str]) -> stim.Tableau:
    X_in,  Z_in  = [stim.PauliString(in_set[i]) for i in range(0,8,2)], [stim.PauliString(in_set[i]) for i in range(1,8,2)]
    X_out, Z_out = [stim.PauliString(out_set[i]) for i in range(0,8,2)], [stim.PauliString(out_set[i]) for i in range(1,8,2)]
    Tin  = stim.Tableau.from_conjugated_generators(
        X_in,
        Z_in
    )
    Tout = stim.Tableau.from_conjugated_generators(
        X_out,
        Z_out
    )
    return Tout @ Tin.inverse()

# ----- пример на твоих множествах -----
IN = [
    "XXZI","XIYZ","YZYI","XYZZ","IXYY","ZYXX","IYIY","ZXZY"
]
OUT = [
    "XIII","YIII","ZXII","ZYII","ZZXI","ZZYI","ZZZX","ZZZY"
]

U = clifford_from_sets(IN, OUT)
print(U.to_circuit())