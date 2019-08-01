#ifndef EXPDATA_H
#define EXPDATA_H

#include <cmath>

struct KAPPAS
{
    double ku;
    double kd;
    double kc;
    double ks;
    double kt;
    double kb;
    double kele;
    double kmuon;
    double ktau;
    double kw;
    double kz;
    double kg;
    double kga;
    double kzga;
};

enum ProcID
{
    iD = 1,
    iU = 2,
    iS = 3,
    iC = 4,
    iB = 5,
    iT = 6,
    iE = 11,
    iM = 13,
    iL = 15,
    iG = 21,
    iA = 22,
    iZ = 23,
    iW = 24,
    iZA = 25,
    // iVBF = 26,
    // iVH = 27,
    iWid = 28
};

double GetKappa(KAPPAS input, int ID);

// The branching ratio and decay width of SM Higgs are obtained from CERN Yellow Report 4 https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
// The branching ratio for those channels not provided in the YR4 are setting to 0
const double bru_SM = 0.0;
const double brd_SM = 0.0;
const double brc_SM = 0.02884;
const double brs_SM = 0.0;
const double brt_SM = 0.0;
const double brb_SM = 0.5809;
const double brele_SM = 0.0;
const double brmuon_SM = 0.0002171;
const double brtau_SM = 0.06256;
const double brW_SM = 0.2152;
const double brZ_SM = 0.02641;
const double brg_SM = 0.0818;
const double brgaga_SM = 0.00227;
const double brZga_SM = 0.001541;
const double H_Width_SM = 0.0041; // in GeV
const double brLeft_SM = 0.0002619;

class SignalStrength
{
public:
    SignalStrength();
    ~SignalStrength(){};
    void SetUpExpData();
    virtual int GetChiSquare(KAPPAS input, double &chisquare);// chisquare is returned in argument, while the true return value is the D.O.F

    static const int NProcMAX = 10;
    int DOF;
    int NProd;
    int NDecay;
    int ProdID[NProcMAX];
    int DecayID[NProcMAX];
    double CentralValue[NProcMAX][NProcMAX];
    double UPError[NProcMAX][NProcMAX];
    double DOError[NProcMAX][NProcMAX];
    bool GoodChannel[NProcMAX][NProcMAX];

};

class STU_EXP
{
public:
    STU_EXP();
    ~STU_EXP(){};
    void SetUpExpData();
    int GetChiSquare(double S, double T, double U, double &chisquare);

    static const int NMAX = 3;
    int DOF;
    int NObs;
    double CentralValue[NMAX];
    double Sigma[NMAX];
    double Rho[NMAX][NMAX];
    double SigmaInv[NMAX][NMAX];
};

class CEPC_STU:public STU_EXP
{
public:
    CEPC_STU();
    ~CEPC_STU(){};
    void SetUpExpData();
};

// class LHC8TeV:public SignalStrength
// {
// // Ref: 1606.02266
// public:
//     LHC8TeV();
//     ~LHC8TeV();
// };

// class ATLAS13TeV:public SignalStrength
// {
// // ATLAS-CONF-2018-031
// // ATLAS-CONF-2016-112
// public:
//     ATLAS13TeV();
//     ~ATLAS13TeV();
// };

// class CMS13TeV:public SignalStrength
// {
// // 1809.10733
// public:
//     CMS13TeV();
//     ~CMS13TeV();
    
// };

// class HLLHC:public SignalStrength
// {
// // JHEP12(2017)153 Table.9
// public:
//     HLLHC();
//     ~HLLHC();
    
// };

// class ILC19:public SignalStrength
// {
// // 1903.01629
// public:
//     ILC19();
//     ~ILC19();
// };

class CEPC:public SignalStrength
{
// 1811.10545
public:
    CEPC();
    ~CEPC(){};
    void SetUpExpData();
};

// class FCCee:public SignalStrength
// {
// // INSPIRE-1713706
// // INSPIRE-1713705
// public:
//     FCCee();
//     ~FCCee();
    
// };
CEPC mu_CEPC;
CEPC_STU STU_CEPC;

void HiggsSignalStrengthSTU_Test(KAPPAS input, double S, double T, double U, double &chi2mu, double &chi2STU, int &DOF, bool &passed);


#endif