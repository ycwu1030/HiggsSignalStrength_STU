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
    iVBF8 = 26,
    iVBF13 = 27,
    iVH13 = 28,
    // iVH = 27,
    iWid = 0
};

const int NEXPmu = 8;
enum EXPmu
{
    muLHC8 = 1,
    muATLAS13 = 2,
    muCMS13 = 4,
    muHLLHC300 = 8,
    muHLLHC3000 = 16,
    muCEPC = 32,
    muILC = 64,
    muFCC = 128 
};

const int NEXPSTU = 4;
enum EXPSTU
{
    STULHC = 1,
    STUCEPC = 2,
    STUILC = 4,
    STUFCC = 8
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
    virtual void SetUpExpData();
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
    virtual void SetUpExpData();
    virtual int GetChiSquare(double S, double T, double U, double &chisquare);

    static const int NMAX = 3;
    int DOF;
    int NObs;
    double CentralValue[NMAX];
    double Sigma[NMAX];
    double Rho[NMAX][NMAX];
    double SigmaInv[NMAX][NMAX];
};

class STU_CEPC:public STU_EXP
{
public:
    STU_CEPC();
    ~STU_CEPC(){};
    virtual void SetUpExpData();
};

class STU_LHC:public STU_EXP
{
public:
    STU_LHC();
    ~STU_LHC(){};
    virtual void SetUpExpData();
    
};

class STU_ILC:public STU_EXP
{
public:
    STU_ILC();
    ~STU_ILC(){};
    virtual void SetUpExpData();
    
};

class STU_FCC:public STU_EXP
{
public:
    STU_FCC();
    ~STU_FCC(){};
    virtual void SetUpExpData();
    
};

class mu_LHC8TeV:public SignalStrength
{
// Ref: 1606.02266
public:
    mu_LHC8TeV();
    ~mu_LHC8TeV(){};
    virtual void SetUpExpData();
};

class mu_ATLAS13TeV:public SignalStrength
{
// ATLAS-CONF-2018-031
// ATLAS-CONF-2016-112
public:
    mu_ATLAS13TeV();
    ~mu_ATLAS13TeV(){};
    virtual void SetUpExpData();
};

class mu_CMS13TeV:public SignalStrength
{
// 1809.10733
public:
    mu_CMS13TeV();
    ~mu_CMS13TeV(){};
    virtual void SetUpExpData();
    
};

class mu_HLLHC300:public SignalStrength
{
// JHEP12(2017)153 Table.9
public:
    mu_HLLHC300();
    ~mu_HLLHC300(){};
    virtual void SetUpExpData();
    virtual int GetChiSquare(KAPPAS input, double &chisquare);
    
private:
    int NChannel[NProcMAX];
    double Contamination[NProcMAX][NProcMAX][NProcMAX];
};

class mu_HLLHC3000:public SignalStrength
{
// JHEP12(2017)153 Table.9
public:
    mu_HLLHC3000();
    ~mu_HLLHC3000(){};
    virtual void SetUpExpData();
    virtual int GetChiSquare(KAPPAS input, double &chisquare);
    
private:
    int NChannel[NProcMAX];
    double Contamination[NProcMAX][NProcMAX][NProcMAX];
};

// class ILC19:public SignalStrength
// {
// // 1903.01629
// public:
//     ILC19();
//     ~ILC19();
// };

class mu_CEPC:public SignalStrength
{
// 1811.10545
public:
    mu_CEPC();
    ~mu_CEPC(){};
    virtual void SetUpExpData();
};

// class FCCee:public SignalStrength
// {
// // INSPIRE-1713706
// // INSPIRE-1713705
// public:
//     FCCee();
//     ~FCCee();
    
// };
// mu_LHC8TeV muExpLHC8;
// mu_ATLAS13TeV muExpATLAS13;
// mu_CMS13TeV muExpCMS13;
// mu_HLLHC300 muExpHLLHC300;
// mu_HLLHC3000 muExpHLLHC3000;
// mu_CEPC muExpCEPC;
// SignalStrength NoExp;

SignalStrength *AllmuExps[NEXPmu] = {new mu_LHC8TeV(), new mu_ATLAS13TeV(), new mu_CMS13TeV(), new mu_HLLHC300(), new mu_HLLHC3000(), new mu_CEPC(), new SignalStrength(), new SignalStrength()};

// STU_LHC STUExpLHC;
// STU_CEPC STUExpCEPC;
// STU_ILC STUExpILC;
// STU_FCC STUExpFCC;

STU_EXP *AllSTUExps[NEXPSTU] = {new STU_LHC(), new STU_CEPC(), new STU_ILC(), new STU_FCC()};

void ChiSquare_Test(double chisquare, int DOF, bool &passed);
void HiggsSignalStrength_Test(int Exps, KAPPAS input, double &chi2mu, int &DOF, bool &passed);
void STU_Test(int Exps, double S, double T, double U, double &chi2STU, int &DOF, bool &passed);
void HiggsSignalStrengthSTU_Test(int muExps, int STUExps, KAPPAS input, double S, double T, double U, double &chi2mu, double &chi2STU, int &DOF, bool &passed);


#endif