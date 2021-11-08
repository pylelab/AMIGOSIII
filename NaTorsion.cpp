const char* docstring=""
"NaTorsion full.pdb [option]\n"
"    calculate torsion angles, bond length and bond angles\n"
"option:\n"
"    1 - (default) torsion angle. missing torsions are marked -360.00.\n"
"        [0] sugar torsion v0: C4'-O4'-C1'-C2'\n"
"        [1] sugar torsion v1: O4'-C1'-C2'-C3'\n"
"        [2] sugar torsion v2: C1'-C2'-C3'-C4'. >0 for C3'-endo; <0 for C2'-endo\n"
"        [3] sugar torsion v3: C2'-C3'-C4'-O4'\n"
"        [4] sugar torsion v4: C3'-C4'-O4'-C1'\n"
"        [5] sugar-backbone torsion v5: C5'-C4'-C3'-C2'\n"
"        [6] sugar-backbone torsion v6: C4'-C3'-C2'-O2'\n"
"        [7] pseudo torsion eta:    C4'[-1]-P-C4'-P[+1]\n"
"        [8] pseudo torsion theta:  P-C4'-P[+1]-C4'[+1]\n"
"        [9] pseudo torsion eta':   C1'[-1]-P-C1'-P[+1]\n"
"        [10]pseudo torsion theta': P-C1'-P[+1]-C1'[+1]\n"
"        [11]backbone torsion alpha:   O3'[-1]-P-O5'-C5'\n"
"        [12]backbone torsion beta:    P-O5'-C5'-C4'\n"
"        [13]backbone torsion gamma:   O5'-C5'-C4'-C3'\n"
"        [14]backbone torsion delta:   C5'-C4'-C3'-O3'\n"
"        [15]backbone torsion epsilon: C4'-C3'-O3'-P[+1]\n"
"        [16]backbone torsion zeta:    C3'-O3'-P[+1]-O5'[+1]\n"
"        [17]improper torsion i1:   C1'-C4'-P-P[+1]\n"
"        [18]sidechain torsion chi: O4'-C1'-N-C\n"
"                                   N=N9,C=C4 for a/g; N=N1,C=C2 for c/t/u\n"
"    2 - bond length. missing bond are marked -1.000\n"
"        [0] pseudo bond: P[-1]-P\n"
"        [1] pseudo bond: O5'[-1]-O5'\n"
"        [2] pseudo bond: C5'[-1]-C5'\n"
"        [3] pseudo bond: C4'[-1]-C4'\n"
"        [4] pseudo bond: C3'[-1]-C3'\n"
"        [5] pseudo bond: O3'[-1]-O3'\n"
"        [6] pseudo bond: C2'[-1]-C2'\n"
"        [7] pseudo bond: C1'[-1]-C1'\n"
"        [8] pseudo bond: O4'[-1]-O4'\n"
"        [9] pseudo bond: P-C4'\n"
"        [10]pseudo bond: C4'-P[+1]\n"
"        [11]pseudo bond: P-C1'\n"
"        [12]pseudo bond: C1'-P[+1]\n"
"        [13]pseudo bond: O3'[-1]-O5'\n"
"        [14]pseudo bond: P-C5'\n"
"        [15]pseudo bond: C5'-C3'\n"
"        [16]pseudo bond: C4'-O3'\n"
"        [17]pseudo bond: C3'-P[+1]\n"
"        [18]pseudo bond: C4'-C2'\n"
"        [19]pseudo bond: O4'-C3'\n"
"        [20]pseudo bond: C4'-C1'\n"
"        [21]pseudo bond: C3'-C1'\n"
"        [22]pseudo bond: O4'-N\n"
"        [23]bond: O3'[-1]-P\n"
"        [24]bond: P-O5'\n"
"        [25]bond: O5'-C5'\n"
"        [26]bond: C5'-C4'\n"
"        [27]bond: C4'-C3'\n"
"        [28]bond: C3'-O3'\n"
"        [29]bond: C3'-C2'\n"
"        [30]bond: C2'-C1'\n"
"        [31]bond: C4'-O4'\n"
"        [32]bond: O4'-C1'\n"
"        [33]bond: C2'-O2'\n"
"        [34]bond: C1'-N\n"
"        [35]bond: N-C\n"
"    4 - bond angle. missing angle are marked -360.00\n"
"        [0] pseudo angle: C4'[-1]-P-C4'\n"
"        [1] pseudo angle: P-C4'-P[+1]\n"
"        [2] pseudo angle: C1'[-1]-P-C1'\n"
"        [3] pseudo angle: P-C1'-P[+1]\n"
"        [4] pseudo angle: P-C4'-C1'\n"
"        [5] angle: O3'[-1]-P-O5'\n"
"        [6] angle: P-O5'-C5'\n"
"        [7] angle: O5'-C5'-C4'\n"
"        [8] angle: C5'-C4'-C3'\n"
"        [9] angle: C4'-C3'-O3'\n"
"        [10]angle: C3'-O3'-P[+1]\n"
"        [11]angle: C4'-C3'-C2'\n"
"        [12]angle: O3'-C3'-C2'\n"
"        [13]angle: C3'-C2'-C1'\n"
"        [14]angle: C3'-C4'-O4'\n"
"        [15]angle: C5'-C4'-O4'\n"
"        [16]angle: C4'-O4'-C1'\n"
"        [17]angle: O4'-C1'-C2'\n"
"        [18]angle: O2'-C2'-C3'\n"
"        [19]angle: O4'-C1'-N\n"
"        [20]angle: C1'-N-C\n"
;

#include <iostream>
#include "PDBParser.hpp"
#include "NaTorsion.hpp"

int main(int argc,char **argv)
{
    if (argc<2)
    {
        cerr<<docstring;
        return 0;
    }
    string filename    =argv[1];
    int option         =(argc>2)?atoi(argv[2]):1;
    bool show_tor      =(option%2==1); option/=2;
    bool show_len      =(option%2==1); option/=2;
    bool show_ang      =(option%2==1);
    
    int atomic_detail  =2;
    int allowX         =1; // only allow ATOM and MSE
    ModelUnit pdb_entry=read_pdb_structure(
        filename.c_str(),atomic_detail,allowX);

    int c,r; // chain index, residue index;
    vector<vector<float> >NaTorMat;
    vector<vector<float> >NaLenMat;
    vector<vector<float> >NaAngMat;
    cout<<"N c resi ";
    if (show_tor)   cout<<"      v0      v1      v2      v3      v4      v5"
        <<"      v6     eta   theta    eta'  theta'   alpha    beta   gamma"
        <<"   delta epsilon    zeta      i1     chi";
    if (show_len)   cout<<"     PmP O5'mO5' C5'mC5' C4'mC4' C3'mC3' O3'mO3'"
        <<" C2'mC2' C1'mC1' O4'mO4'    PC4'   C4'Pp    PC1'   C1'Pp O3'mO5'"
        <<"    PC5'  C5'C3'  C4'O3'   C3'Pp  C4'C2'  O4'C3'  C4'C1'  C3'C1'"
        <<"    O4'N  O3'mP    PO5'   O5'C5'  C5'C4'  C4'C3'  C3'O3'  C3'C2'"
        <<"  C2'C1'  C4'O4'  O4'C1'  C2'O2'    C1'N      NC";
    if (show_ang)       cout<<"  C4'mPC4'    PC4'Pp  C1'mPC1'    PC1'Pp   PC4'C1'"
        <<"  O3'mPO5'   PO5'C5' O5'C5'C4' C5'C4'C3' C4'C3'O3'  C3'O3'Pp C4'C3'C2'"
        <<" O3'C3'C2' C3'C2'C1' C3'C4'O4' C5'C4'O4' C4'O4'C1' O4'C1'C2' O2'C2'C3'"
        <<"   O4'C1'N     C1'NC";
    cout<<endl;
    for (c=0;c<pdb_entry.chains.size();c++)
    {
        NaTorsion(pdb_entry.chains[c],NaTorMat,NaLenMat,NaAngMat,
            show_tor,show_len,show_ang);
        for (r=0;r<pdb_entry.chains[c].residues.size();r++)
        {
            cout<<char(tolower(pdb_entry.chains[c].residues[r].resn[2]))<<' '
                <<pdb_entry.chains[c].chainID<<' '
                <<setw(4)<<pdb_entry.chains[c].residues[r].resi
                <<pdb_entry.chains[c].residues[r].icode;
            if (show_tor)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(2)
                    <<' '<<setw(7)<<NaTorMat[r][0]
                    <<' '<<setw(7)<<NaTorMat[r][1]
                    <<' '<<setw(7)<<NaTorMat[r][2]
                    <<' '<<setw(7)<<NaTorMat[r][3]
                    <<' '<<setw(7)<<NaTorMat[r][4]
                    <<' '<<setw(7)<<NaTorMat[r][5]
                    <<' '<<setw(7)<<NaTorMat[r][6]
                    <<' '<<setw(7)<<NaTorMat[r][7]
                    <<' '<<setw(7)<<NaTorMat[r][8]
                    <<' '<<setw(7)<<NaTorMat[r][9]
                    <<' '<<setw(7)<<NaTorMat[r][10]
                    <<' '<<setw(7)<<NaTorMat[r][11]
                    <<' '<<setw(7)<<NaTorMat[r][12]
                    <<' '<<setw(7)<<NaTorMat[r][13]
                    <<' '<<setw(7)<<NaTorMat[r][14]
                    <<' '<<setw(7)<<NaTorMat[r][15]
                    <<' '<<setw(7)<<NaTorMat[r][16]
                    <<' '<<setw(7)<<NaTorMat[r][17]
                    <<' '<<setw(7)<<NaTorMat[r][18];
                NaTorMat[r].clear();
            }
            if (show_len)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(3)
                    <<' '<<setw(7)<<NaLenMat[r][0]
                    <<' '<<setw(7)<<NaLenMat[r][1]
                    <<' '<<setw(7)<<NaLenMat[r][2]
                    <<' '<<setw(7)<<NaLenMat[r][3]
                    <<' '<<setw(7)<<NaLenMat[r][4]
                    <<' '<<setw(7)<<NaLenMat[r][5]
                    <<' '<<setw(7)<<NaLenMat[r][6]
                    <<' '<<setw(7)<<NaLenMat[r][7]
                    <<' '<<setw(7)<<NaLenMat[r][8]
                    <<' '<<setw(7)<<NaLenMat[r][9]
                    <<' '<<setw(7)<<NaLenMat[r][10]
                    <<' '<<setw(7)<<NaLenMat[r][11]
                    <<' '<<setw(7)<<NaLenMat[r][12]
                    <<' '<<setw(7)<<NaLenMat[r][13]
                    <<' '<<setw(7)<<NaLenMat[r][14]
                    <<' '<<setw(7)<<NaLenMat[r][15]
                    <<' '<<setw(7)<<NaLenMat[r][16]
                    <<' '<<setw(7)<<NaLenMat[r][17]
                    <<' '<<setw(7)<<NaLenMat[r][18]
                    <<' '<<setw(7)<<NaLenMat[r][19]
                    <<' '<<setw(7)<<NaLenMat[r][20]
                    <<' '<<setw(7)<<NaLenMat[r][21]
                    <<' '<<setw(7)<<NaLenMat[r][22]
                    <<' '<<setw(7)<<NaLenMat[r][23]
                    <<' '<<setw(7)<<NaLenMat[r][24]
                    <<' '<<setw(7)<<NaLenMat[r][25]
                    <<' '<<setw(7)<<NaLenMat[r][26]
                    <<' '<<setw(7)<<NaLenMat[r][27]
                    <<' '<<setw(7)<<NaLenMat[r][28]
                    <<' '<<setw(7)<<NaLenMat[r][29]
                    <<' '<<setw(7)<<NaLenMat[r][30]
                    <<' '<<setw(7)<<NaLenMat[r][31]
                    <<' '<<setw(7)<<NaLenMat[r][32]
                    <<' '<<setw(7)<<NaLenMat[r][33]
                    <<' '<<setw(7)<<NaLenMat[r][34]
                    <<' '<<setw(7)<<NaLenMat[r][35];
                    NaLenMat[r].clear();
            }
            if (show_ang)
            {
                cout<<setiosflags(ios::fixed)<<setprecision(2)
                    <<' '<<setw(9)<<NaAngMat[r][0]
                    <<' '<<setw(9)<<NaAngMat[r][1]
                    <<' '<<setw(9)<<NaAngMat[r][2]
                    <<' '<<setw(9)<<NaAngMat[r][3]
                    <<' '<<setw(9)<<NaAngMat[r][4]
                    <<' '<<setw(9)<<NaAngMat[r][5]
                    <<' '<<setw(9)<<NaAngMat[r][6]
                    <<' '<<setw(9)<<NaAngMat[r][7]
                    <<' '<<setw(9)<<NaAngMat[r][8]
                    <<' '<<setw(9)<<NaAngMat[r][9]
                    <<' '<<setw(9)<<NaAngMat[r][10]
                    <<' '<<setw(9)<<NaAngMat[r][11]
                    <<' '<<setw(9)<<NaAngMat[r][12]
                    <<' '<<setw(9)<<NaAngMat[r][13]
                    <<' '<<setw(9)<<NaAngMat[r][14]
                    <<' '<<setw(9)<<NaAngMat[r][15]
                    <<' '<<setw(9)<<NaAngMat[r][16]
                    <<' '<<setw(9)<<NaAngMat[r][17]
                    <<' '<<setw(9)<<NaAngMat[r][18]
                    <<' '<<setw(9)<<NaAngMat[r][19]
                    <<' '<<setw(9)<<NaAngMat[r][20];
                    NaAngMat[r].clear();
            }
            cout<<endl;
        }
        if (show_tor) NaTorMat.clear();
        if (show_len) NaLenMat.clear();
        if (show_ang) NaAngMat.clear();
    }
    return 0;
}
