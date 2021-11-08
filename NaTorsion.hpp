/* calculate backbone torsion angles (omega, psi, phi, kappa, alpha) */
#ifndef NaTorsion_HPP
#define NaTorsion_HPP 1
#include <cstring>
#include "PDBParser.hpp"
#include "GeometryTools.hpp"

using namespace std;


/* calculate pucker (C1'C2'C3'C4'), pseudo torsion (eta, theta, eta', theta')
 * backbone torsion (alpha, beta, gamma, delta , epsilon, zeta) and
 * sidechain torsion (chi) from backbone.pdb */
void NaTorsion(ChainUnit& chain, vector<vector<float> >&NaTorMat,
    vector<vector<float> >&NaLenMat, vector<vector<float> >&NaAngMat, 
    const bool show_tor, const bool show_len, const bool show_ang)
{
    int L=chain.residues.size();
    // default torsion angles: omega=phi=psi=360
    vector<float> tmp_tor(19,-360.);
    vector<float> tmp_len(36,-1.);
    vector<float> tmp_ang(21,-360.);
    if (show_tor) NaTorMat.assign(L,tmp_tor);
    if (show_len) NaLenMat.assign(L,tmp_len);
    if (show_ang) NaAngMat.assign(L,tmp_ang);

    // coordinates of previous residue
    vector<float> prev_P(3,0);  bool has_prev_P=false;
    vector<float> prev_O5(3,0); bool has_prev_O5=false;
    vector<float> prev_C5(3,0); bool has_prev_C5=false;
    vector<float> prev_C4(3,0); bool has_prev_C4=false;
    vector<float> prev_C3(3,0); bool has_prev_C3=false;
    vector<float> prev_C2(3,0); bool has_prev_C2=false;
    vector<float> prev_C1(3,0); bool has_prev_C1=false;
    vector<float> prev_O4(3,0); bool has_prev_O4=false;
    vector<float> prev_O3(3,0); bool has_prev_O3=false;
    // coordinates of current residue
    vector<float> P(3,0.);      bool has_P=false;
    vector<float> O5(3,0.);     bool has_O5=false;
    vector<float> C5(3,0.);     bool has_C5=false;
    vector<float> C4(3,0.);     bool has_C4=false;
    vector<float> C3(3,0.);     bool has_C3=false;
    vector<float> C2(3,0.);     bool has_C2=false;
    vector<float> C1(3,0.);     bool has_C1=false;
    vector<float> O4(3,0.);     bool has_O4=false;
    vector<float> O3(3,0.);     bool has_O3=false;
    vector<float> O2(3,0.);     bool has_O2=false;
    vector<float> Nx(3,0.);     bool has_Nx=false;
    vector<float> Cx(3,0.);     bool has_Cx=false;
    // coordinates of next residue
    vector<float> next_P(3,0);  bool has_next_P=false;
    vector<float> next_O5(3,0); bool has_next_O5=false;
    vector<float> next_C4(3,0); bool has_next_C4=false;
    vector<float> next_C1(3,0); bool has_next_C1=false;
    char base=' ';

    int r,a; // residue index, atom index
    
    for (r=0;r<L;r++) 
    {
        // whether the required atom exists
        has_prev_P=false;
        has_prev_O5=false;
        has_prev_C5=false;
        has_prev_C4=false;
        has_prev_C3=false;
        has_prev_C2=false;
        has_prev_C1=false;
        has_prev_O4=false;
        has_prev_O3=false;
        has_P=false;
        has_O5=false;
        has_C5=false;
        has_C4=false;
        has_C3=false;
        has_C2=false;
        has_C1=false;
        has_O4=false;
        has_O3=false;
        has_O2=false;
        has_Nx=false;
        has_Cx=false;
        has_next_P=false;
        has_next_O5=false;
        has_next_C4=false;
        has_next_C1=false;
        base=tolower(chain.residues[r].resn[2]);
        
        if (r>0) // find previous residue atoms
        {
            for (a=0;a<chain.residues[r-1].atoms.size();a++)
            {
                if (chain.residues[r-1].atoms[a].name==" P  ")
                {
                    has_prev_P=true;
                    prev_P=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" O5'")
                {
                    has_prev_O5=true;
                    prev_O5=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C5'")
                {
                    has_prev_C5=true;
                    prev_C5=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C4'")
                {
                    has_prev_C4=true;
                    prev_C4=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C3'")
                {
                    has_prev_C3=true;
                    prev_C3=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C2'")
                {
                    has_prev_C2=true;
                    prev_C2=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" C1'")
                {
                    has_prev_C1=true;
                    prev_C1=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" O4'")
                {
                    has_prev_O4=true;
                    prev_O4=chain.residues[r-1].atoms[a].xyz;
                }
                else if (chain.residues[r-1].atoms[a].name==" O3'")
                {
                    has_prev_O3=true;
                    prev_O3=chain.residues[r-1].atoms[a].xyz;
                }
            }
        }

        // find current residue atoms
        for (a=0;a<chain.residues[r].atoms.size();a++)
        {
            if (chain.residues[r].atoms[a].name==" P  ")
            {
                has_P=true;
                P=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O5'")
            {
                has_O5=true;
                O5=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C5'")
            {
                has_C5=true;
                C5=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C4'")
            {
                has_C4=true;
                C4=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C3'")
            {
                has_C3=true;
                C3=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C2'")
            {
                has_C2=true;
                C2=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" C1'")
            {
                has_C1=true;
                C1=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O4'")
            {
                has_O4=true;
                O4=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O3'")
            {
                has_O3=true;
                O3=chain.residues[r].atoms[a].xyz;
            }
            else if (chain.residues[r].atoms[a].name==" O2'")
            {
                has_O2=true;
                O2=chain.residues[r].atoms[a].xyz;
            }
            else if((chain.residues[r].atoms[a].name==" N9 " &&
                    (base=='a' || base=='g')) ||
                    (chain.residues[r].atoms[a].name==" N1 " && 
                    (base=='c' || base=='t' || base=='u')))
            {
                has_Nx=true;
                Nx=chain.residues[r].atoms[a].xyz;
            }
            else if((chain.residues[r].atoms[a].name==" C4 " &&
                    (base=='a' || base=='g')) ||
                    (chain.residues[r].atoms[a].name==" C2 " && 
                    (base=='c' || base=='t' || base=='u')))
            {
                has_Cx=true;
                Cx=chain.residues[r].atoms[a].xyz;
            }
        }

        if (r+1<L) // find next residue atoms
        {
            for (a=0;a<chain.residues[r+1].atoms.size();a++)
            {
                if (chain.residues[r+1].atoms[a].name==" P  ")
                {
                    has_next_P=true;
                    next_P=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" O5'")
                {
                    has_next_O5=true;
                    next_O5=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" C4'")
                {
                    has_next_C4=true;
                    next_C4=chain.residues[r+1].atoms[a].xyz;
                }
                else if (chain.residues[r+1].atoms[a].name==" C1'")
                {
                    has_next_C1=true;
                    next_C1=chain.residues[r+1].atoms[a].xyz;
                }
            }
        }
 
        if (show_tor)
        {
            if (has_C4 && has_O4 && has_C1 && has_C2)         NaTorMat[r][0]=rad2deg(Points2Dihedral(C4,O4,C1,C2));          // v0
            if (has_O4 && has_C1 && has_C2 && has_C3)         NaTorMat[r][1]=rad2deg(Points2Dihedral(O4,C1,C2,C3));          // v1
            if (has_C1 && has_C2 && has_C3 && has_C4)         NaTorMat[r][2]=rad2deg(Points2Dihedral(C1,C2,C3,C4));          // v2 pucker
            if (has_C2 && has_C3 && has_C4 && has_O4)         NaTorMat[r][3]=rad2deg(Points2Dihedral(C2,C3,C4,O4));          // v3
            if (has_C3 && has_C4 && has_O4 && has_C1)         NaTorMat[r][4]=rad2deg(Points2Dihedral(C3,C4,O4,C1));          // v4
            if (has_C5 && has_C4 && has_C3 && has_C2)         NaTorMat[r][5]=rad2deg(Points2Dihedral(C5,C4,C3,C2));          // v5
            if (has_C4 && has_C3 && has_C2 && has_O2)         NaTorMat[r][6]=rad2deg(Points2Dihedral(C4,C3,C2,O2));          // v6
            if (has_prev_C4 && has_P && has_C4 && has_next_P) NaTorMat[r][7]=rad2deg(Points2Dihedral(prev_C4,P,C4,next_P));  // eta
            if (has_P && has_C4 && has_next_P && has_next_C4) NaTorMat[r][8]=rad2deg(Points2Dihedral(P,C4,next_P,next_C4));  // theta
            if (has_prev_C1 && has_P && has_C1 && has_next_P) NaTorMat[r][9]=rad2deg(Points2Dihedral(prev_C1,P,C1,next_P));  // eta'
            if (has_P && has_C1 && has_next_P && has_next_C1) NaTorMat[r][10]=rad2deg(Points2Dihedral(P,C1,next_P,next_C1));  // theta'
            if (has_prev_O3 && has_P && has_O5 && has_C5)     NaTorMat[r][11]=rad2deg(Points2Dihedral(prev_O3,P,O5,C5));     // alpha
            if (has_P  && has_O5 && has_C5 && has_C4)         NaTorMat[r][12]=rad2deg(Points2Dihedral(P,O5,C5,C4));          // beta
            if (has_O5 && has_C5 && has_C4 && has_C3)         NaTorMat[r][13]=rad2deg(Points2Dihedral(O5,C5,C4,C3));         // gamma
            if (has_C5 && has_C4 && has_C3 && has_O3)         NaTorMat[r][14]=rad2deg(Points2Dihedral(C5,C4,C3,O3));         // delta
            if (has_C4 && has_C3 && has_O3 && has_next_P)     NaTorMat[r][15]=rad2deg(Points2Dihedral(C4,C3,O3,next_P));     // epsilon
            if (has_C3 && has_O3 && has_next_P && has_next_O5)NaTorMat[r][16]=rad2deg(Points2Dihedral(C3,O3,next_P,next_O5));// zeta
            if (has_C1 && has_C4 && has_P && has_next_P)      NaTorMat[r][17]=rad2deg(Points2Dihedral(C1,C4,P,next_P));      // i1
            if (has_O4 && has_C1 && has_Nx && has_Cx)         NaTorMat[r][18]=rad2deg(Points2Dihedral(O4,C1,Nx,Cx));         // chi
        }
        if (show_len)
        {
            if (has_prev_P  && has_P)  NaLenMat[r][0] =Points2Distance(prev_P,P);   // P[-1]-P
            if (has_prev_O5 && has_O5) NaLenMat[r][1] =Points2Distance(prev_O5,O5); // O5'[-1]-O5'
            if (has_prev_C5 && has_C5) NaLenMat[r][2] =Points2Distance(prev_C5,C5); // C5'[-1]-C5'
            if (has_prev_C4 && has_C4) NaLenMat[r][3] =Points2Distance(prev_C4,C4); // C4'[-1]-C4'
            if (has_prev_C3 && has_C3) NaLenMat[r][4] =Points2Distance(prev_C3,C3); // C3'[-1]-C3'
            if (has_prev_O3 && has_O3) NaLenMat[r][5] =Points2Distance(prev_O3,O3); // O3'[-1]-O3'
            if (has_prev_C2 && has_C2) NaLenMat[r][6] =Points2Distance(prev_C2,C2); // C2'[-1]-C2'
            if (has_prev_C1 && has_C1) NaLenMat[r][7] =Points2Distance(prev_C1,C1); // C1'[-1]-C1'
            if (has_prev_O4 && has_O4) NaLenMat[r][8] =Points2Distance(prev_O4,O4); // O4'[-1]-O4'
            if (has_P       && has_C4) NaLenMat[r][9] =Points2Distance(P,C4);       // P-C4'
            if (has_C4  && has_next_P) NaLenMat[r][10]=Points2Distance(C4,next_P);  // C4'-P[+1]
            if (has_P       && has_C1) NaLenMat[r][11]=Points2Distance(P,C1);       // P-C1'
            if (has_C1  && has_next_P) NaLenMat[r][12]=Points2Distance(C1,next_P);  // C1'-P[+1]
            if (has_prev_O3 && has_O5) NaLenMat[r][13]=Points2Distance(prev_O3,O5); // O3'[-1]-O5'
            if (has_P       && has_C5) NaLenMat[r][14]=Points2Distance(P,C5);       // P-C5'
            if (has_C5      && has_C3) NaLenMat[r][15]=Points2Distance(C5,C3);      // C5'-C3'
            if (has_C4      && has_O3) NaLenMat[r][16]=Points2Distance(C4,O3);      // C4'-O3'
            if (has_C3  && has_next_P) NaLenMat[r][17]=Points2Distance(C3,next_P);  // C3'-P[+1]
            if (has_C4      && has_C2) NaLenMat[r][18]=Points2Distance(C4,C2);      // C4'-C2'
            if (has_O4      && has_C3) NaLenMat[r][19]=Points2Distance(O4,C3);      // O4'-C3'
            if (has_C4      && has_C1) NaLenMat[r][20]=Points2Distance(C4,C1);      // C4'-C1'
            if (has_C3      && has_C1) NaLenMat[r][21]=Points2Distance(C3,C1);      // C3'-C1'
            if (has_O4      && has_Nx) NaLenMat[r][22]=Points2Distance(O4,Nx);      // O4'-N
            if (has_prev_O3 && has_P)  NaLenMat[r][23]=Points2Distance(prev_O3,P);  // O3'[-1]-P
            if (has_P       && has_O5) NaLenMat[r][24]=Points2Distance(P,O5);       // P-O5'
            if (has_O5      && has_C5) NaLenMat[r][25]=Points2Distance(O5,C5);      // O5'-C5'
            if (has_C5      && has_C4) NaLenMat[r][26]=Points2Distance(C5,C4);      // C5'-C4'
            if (has_C4      && has_C3) NaLenMat[r][27]=Points2Distance(C4,C3);      // C4'-C3'
            if (has_C3      && has_O3) NaLenMat[r][28]=Points2Distance(C3,O3);      // C3'-O3'
            if (has_C3      && has_C2) NaLenMat[r][29]=Points2Distance(C3,C2);      // C3'-C2'
            if (has_C2      && has_C1) NaLenMat[r][30]=Points2Distance(C2,C1);      // C2'-C1'
            if (has_C4      && has_O4) NaLenMat[r][31]=Points2Distance(C4,O4);      // C4'-O4'
            if (has_O4      && has_C1) NaLenMat[r][32]=Points2Distance(O4,C1);      // O4'-O1'
            if (has_C2      && has_O2) NaLenMat[r][33]=Points2Distance(C2,O2);      // C2'-O2'
            if (has_C1      && has_Nx) NaLenMat[r][34]=Points2Distance(C1,Nx);      // C1'-Nx
            if (has_Nx      && has_Cx) NaLenMat[r][35]=Points2Distance(Nx,Cx);      // N-C
        }
        if (show_ang)
        {
            if (has_prev_C4 && has_P && has_C4) NaAngMat[r][0] =rad2deg(Points2Angle(prev_C4,P,C4)); // C4'[-1]-P-C4'
            if (has_P  && has_C4 && has_next_P) NaAngMat[r][1] =rad2deg(Points2Angle(P,C4,next_P));  // P-C4'-P[+1]
            if (has_prev_C1 && has_P && has_C1) NaAngMat[r][2] =rad2deg(Points2Angle(prev_C1,P,C1)); // C1'[-1]-P-C1'
            if (has_P  && has_C1 && has_next_P) NaAngMat[r][3] =rad2deg(Points2Angle(P,C1,next_P));  // P-C1'-P[+1]
            if (has_P  && has_C4 && has_C1)     NaAngMat[r][4] =rad2deg(Points2Angle(P,C4,C1));      // P-C4'-C1'
            if (has_prev_O3 && has_P && has_O5) NaAngMat[r][5] =rad2deg(Points2Angle(prev_O3,P,O5)); // O3'[-1]-P-O5'
            if (has_P  && has_O5 && has_C5)     NaAngMat[r][6] =rad2deg(Points2Angle(P,O5,C5));      // P-O5'-C5'
            if (has_O5 && has_C5 && has_C4)     NaAngMat[r][7] =rad2deg(Points2Angle(O5,C5,C4));     // O5'-C5'-C4'
            if (has_C5 && has_C4 && has_C3)     NaAngMat[r][8] =rad2deg(Points2Angle(C5,C4,C3));     // C5'-C4'-C3'
            if (has_C4 && has_C3 && has_O3)     NaAngMat[r][9] =rad2deg(Points2Angle(C4,C3,O3));     // C4'-C3'-O3'
            if (has_C3 && has_O3 && has_next_P) NaAngMat[r][10]=rad2deg(Points2Angle(C3,O3,next_P)); // C3'-O3'-P[+1]
            if (has_C4 && has_C3 && has_C2)     NaAngMat[r][11]=rad2deg(Points2Angle(C4,C3,C2));     // C4'-C3'-C2'
            if (has_O3 && has_C3 && has_C2)     NaAngMat[r][12]=rad2deg(Points2Angle(O3,C3,C2));     // O3'-C3'-C2'
            if (has_C3 && has_C2 && has_C1)     NaAngMat[r][13]=rad2deg(Points2Angle(C3,C2,C1));     // C3'-C2'-C1'
            if (has_C3 && has_C4 && has_O4)     NaAngMat[r][14]=rad2deg(Points2Angle(C3,C4,O4));     // C3'-C4'-O4'
            if (has_C5 && has_C4 && has_O4)     NaAngMat[r][15]=rad2deg(Points2Angle(C5,C4,O4));     // C5'-C4'-O4'
            if (has_C4 && has_O4 && has_C1)     NaAngMat[r][16]=rad2deg(Points2Angle(C4,O4,C1));     // C4'-O4'-C1'
            if (has_O4 && has_C1 && has_C2)     NaAngMat[r][17]=rad2deg(Points2Angle(O4,C1,C2));     // O4'-C1'-C2'
            if (has_O2 && has_C2 && has_C3)     NaAngMat[r][18]=rad2deg(Points2Angle(O2,C2,C3));     // O2'-C2'-C3'
            if (has_O4 && has_C1 && has_Nx)     NaAngMat[r][19]=rad2deg(Points2Angle(O4,C1,Nx));     // O4'-C1'-Nx
            if (has_C1 && has_Nx && has_Cx)     NaAngMat[r][20]=rad2deg(Points2Angle(C1,Nx,Cx));     // C1'-Nx-Cx
        }
    }
    

    tmp_tor.clear();
    tmp_len.clear();
    tmp_ang.clear();
    prev_C4.clear();
    prev_C1.clear();
    prev_O3.clear();
    P.clear();
    O5.clear();
    C5.clear();
    C4.clear();
    C3.clear();
    C2.clear();
    C1.clear();
    O4.clear();
    O3.clear();
    Nx.clear();
    Cx.clear();
    next_P.clear(); 
    next_O5.clear();
    next_C4.clear();
    next_C1.clear();
    return;
}

#endif
