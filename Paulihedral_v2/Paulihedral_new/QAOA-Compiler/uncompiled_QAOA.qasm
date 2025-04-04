OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
creg c[20];
h q[16];
h q[18];
h q[17];
h q[15];
h q[14];
h q[13];
h q[12];
h q[19];
h q[11];
h q[10];
h q[9];
h q[8];
h q[7];
h q[6];
h q[5];
h q[4];
h q[3];
h q[2];
h q[1];
h q[0];
cx q[16],q[18];
rz(g1_16_18) q[18];
cx q[16],q[18];
cx q[16],q[17];
rz(g1_16_17) q[17];
cx q[16],q[17];
cx q[16],q[15];
rz(g1_16_15) q[15];
cx q[16],q[15];
cx q[16],q[14];
rz(g1_16_14) q[14];
cx q[16],q[14];
cx q[16],q[9];
rz(g1_16_9) q[9];
cx q[16],q[9];
cx q[16],q[8];
rz(g1_16_8) q[8];
cx q[16],q[8];
cx q[16],q[7];
rz(g1_16_7) q[7];
cx q[16],q[7];
cx q[16],q[4];
rz(g1_16_4) q[4];
cx q[16],q[4];
cx q[16],q[3];
rz(g1_16_3) q[3];
cx q[16],q[3];
cx q[16],q[1];
rz(g1_16_1) q[1];
cx q[16],q[1];
cx q[16],q[0];
rz(g1_16_0) q[0];
cx q[16],q[0];
cx q[18],q[14];
rz(g1_18_14) q[14];
cx q[18],q[14];
cx q[18],q[13];
rz(g1_18_13) q[13];
cx q[18],q[13];
cx q[18],q[11];
rz(g1_18_11) q[11];
cx q[18],q[11];
cx q[18],q[5];
rz(g1_18_5) q[5];
cx q[18],q[5];
cx q[18],q[3];
rz(g1_18_3) q[3];
cx q[18],q[3];
cx q[17],q[14];
rz(g1_17_14) q[14];
cx q[17],q[14];
cx q[17],q[13];
rz(g1_17_13) q[13];
cx q[17],q[13];
cx q[17],q[11];
rz(g1_17_11) q[11];
cx q[17],q[11];
cx q[17],q[7];
rz(g1_17_7) q[7];
cx q[17],q[7];
cx q[17],q[5];
rz(g1_17_5) q[5];
cx q[17],q[5];
cx q[17],q[4];
rz(g1_17_4) q[4];
cx q[17],q[4];
cx q[17],q[2];
rz(g1_17_2) q[2];
cx q[17],q[2];
cx q[15],q[14];
rz(g1_15_14) q[14];
cx q[15],q[14];
cx q[15],q[11];
rz(g1_15_11) q[11];
cx q[15],q[11];
cx q[15],q[3];
rz(g1_15_3) q[3];
cx q[15],q[3];
cx q[15],q[1];
rz(g1_15_1) q[1];
cx q[15],q[1];
cx q[15],q[0];
rz(g1_15_0) q[0];
cx q[15],q[0];
cx q[14],q[12];
rz(g1_14_12) q[12];
cx q[14],q[12];
cx q[14],q[11];
rz(g1_14_11) q[11];
cx q[14],q[11];
cx q[14],q[9];
rz(g1_14_9) q[9];
cx q[14],q[9];
cx q[14],q[8];
rz(g1_14_8) q[8];
cx q[14],q[8];
cx q[14],q[5];
rz(g1_14_5) q[5];
cx q[14],q[5];
cx q[14],q[4];
rz(g1_14_4) q[4];
cx q[14],q[4];
cx q[14],q[1];
rz(g1_14_1) q[1];
cx q[14],q[1];
cx q[14],q[0];
rz(g1_14_0) q[0];
cx q[14],q[0];
cx q[13],q[11];
rz(g1_13_11) q[11];
cx q[13],q[11];
cx q[13],q[10];
rz(g1_13_10) q[10];
cx q[13],q[10];
cx q[13],q[7];
rz(g1_13_7) q[7];
cx q[13],q[7];
cx q[13],q[4];
rz(g1_13_4) q[4];
cx q[13],q[4];
cx q[13],q[3];
rz(g1_13_3) q[3];
cx q[13],q[3];
cx q[13],q[1];
rz(g1_13_1) q[1];
cx q[13],q[1];
cx q[13],q[0];
rz(g1_13_0) q[0];
cx q[13],q[0];
cx q[12],q[19];
rz(g1_12_19) q[19];
cx q[12],q[19];
cx q[12],q[11];
rz(g1_12_11) q[11];
cx q[12],q[11];
cx q[12],q[10];
rz(g1_12_10) q[10];
cx q[12],q[10];
cx q[12],q[8];
rz(g1_12_8) q[8];
cx q[12],q[8];
cx q[12],q[3];
rz(g1_12_3) q[3];
cx q[12],q[3];
cx q[12],q[2];
rz(g1_12_2) q[2];
cx q[12],q[2];
cx q[19],q[11];
rz(g1_19_11) q[11];
cx q[19],q[11];
cx q[19],q[6];
rz(g1_19_6) q[6];
cx q[19],q[6];
cx q[19],q[5];
rz(g1_19_5) q[5];
cx q[19],q[5];
cx q[19],q[2];
rz(g1_19_2) q[2];
cx q[19],q[2];
cx q[19],q[1];
rz(g1_19_1) q[1];
cx q[19],q[1];
cx q[19],q[0];
rz(g1_19_0) q[0];
cx q[19],q[0];
cx q[11],q[10];
rz(g1_11_10) q[10];
cx q[11],q[10];
cx q[11],q[7];
rz(g1_11_7) q[7];
cx q[11],q[7];
cx q[11],q[5];
rz(g1_11_5) q[5];
cx q[11],q[5];
cx q[11],q[3];
rz(g1_11_3) q[3];
cx q[11],q[3];
cx q[11],q[1];
rz(g1_11_1) q[1];
cx q[11],q[1];
cx q[11],q[0];
rz(g1_11_0) q[0];
cx q[11],q[0];
cx q[10],q[8];
rz(g1_10_8) q[8];
cx q[10],q[8];
cx q[10],q[5];
rz(g1_10_5) q[5];
cx q[10],q[5];
cx q[10],q[3];
rz(g1_10_3) q[3];
cx q[10],q[3];
cx q[10],q[2];
rz(g1_10_2) q[2];
cx q[10],q[2];
cx q[10],q[0];
rz(g1_10_0) q[0];
cx q[10],q[0];
cx q[9],q[5];
rz(g1_9_5) q[5];
cx q[9],q[5];
cx q[8],q[7];
rz(g1_8_7) q[7];
cx q[8],q[7];
cx q[8],q[6];
rz(g1_8_6) q[6];
cx q[8],q[6];
cx q[8],q[4];
rz(g1_8_4) q[4];
cx q[8],q[4];
cx q[8],q[1];
rz(g1_8_1) q[1];
cx q[8],q[1];
cx q[7],q[5];
rz(g1_7_5) q[5];
cx q[7],q[5];
cx q[7],q[4];
rz(g1_7_4) q[4];
cx q[7],q[4];
cx q[7],q[3];
rz(g1_7_3) q[3];
cx q[7],q[3];
cx q[7],q[1];
rz(g1_7_1) q[1];
cx q[7],q[1];
cx q[6],q[5];
rz(g1_6_5) q[5];
cx q[6],q[5];
cx q[6],q[4];
rz(g1_6_4) q[4];
cx q[6],q[4];
cx q[6],q[1];
rz(g1_6_1) q[1];
cx q[6],q[1];
cx q[6],q[0];
rz(g1_6_0) q[0];
cx q[6],q[0];
cx q[5],q[1];
rz(g1_5_1) q[1];
cx q[5],q[1];
cx q[5],q[0];
rz(g1_5_0) q[0];
cx q[5],q[0];
cx q[4],q[3];
rz(g1_4_3) q[3];
cx q[4],q[3];
cx q[3],q[2];
rz(g1_3_2) q[2];
cx q[3],q[2];
cx q[2],q[0];
rz(g1_2_0) q[0];
cx q[2],q[0];
rx(b1) q[16];
rx(b1) q[18];
rx(b1) q[17];
rx(b1) q[15];
rx(b1) q[14];
rx(b1) q[13];
rx(b1) q[12];
rx(b1) q[19];
rx(b1) q[11];
rx(b1) q[10];
rx(b1) q[9];
rx(b1) q[8];
rx(b1) q[7];
rx(b1) q[6];
rx(b1) q[5];
rx(b1) q[4];
rx(b1) q[3];
rx(b1) q[2];
rx(b1) q[1];
rx(b1) q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15],q[16],q[17],q[18],q[19];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
measure q[17] -> c[17];
measure q[18] -> c[18];
measure q[19] -> c[19];
