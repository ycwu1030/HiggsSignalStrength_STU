#include "HSS_STU_ExpData.h"
#include <iostream>

SignalStrength::SignalStrength()
{
    SetUpExpData();
}
void SignalStrength::SetUpExpData()
{
    NProd = 0;
    NDecay = 0;
}
int SignalStrength::GetChiSquare(KAPPAS input, double &chisquare)
{
    chisquare = 0.0;
    DOF = 0;
    double kwid = GetKappa(input,iWid);
    double kprod;
    double kdecay;
    double mu;
    for (int ipro = 0; ipro < NProd; ++ipro)
    {
        for (int idec = 0; idec < NDecay; ++idec)
        {
            if (!GoodChannel[ipro][idec]) continue;
            DOF+=1;
            kprod = GetKappa(input,ProdID[ipro]);
            kdecay = GetKappa(input,DecayID[idec]);
            mu = pow(kprod*kdecay,2)/kwid;
            if (mu >= CentralValue[ipro][idec])
            {
                chisquare+=pow((mu - CentralValue[ipro][idec])/UPError[ipro][idec],2);
            }
            else
            {
                chisquare+=pow((mu - CentralValue[ipro][idec])/DOError[ipro][idec],2);
            }
            // std::cout<<"ipro "<<ipro<<" idec "<<idec<<" chisquare "<<chisquare<<std::endl;
        }
    }
    return DOF;
}

mu_CEPC::mu_CEPC()
{
    SetUpExpData();
}

void mu_CEPC::SetUpExpData()
{
    NProd = 1;
    NDecay = 8;
    ProdID[0] = iZ;
    DecayID[0] = iB;
    DecayID[1] = iC;
    DecayID[2] = iG;
    DecayID[3] = iW;
    DecayID[4] = iL;
    DecayID[5] = iZ;
    DecayID[6] = iA;
    DecayID[7] = iM;
    UPError[0][0] = 0.0027;
    UPError[0][1] = 0.033;
    UPError[0][2] = 0.013;
    UPError[0][3] = 0.01;
    UPError[0][4] = 0.008;
    UPError[0][5] = 0.051;
    UPError[0][6] = 0.068;
    UPError[0][7] = 0.17;
    for (int ipro = 0; ipro < NProd; ++ipro)
    {
        for (int idec = 0; idec < NDecay; ++idec)
        {
            CentralValue[ipro][idec] = 1.0;
            DOError[ipro][idec] = UPError[ipro][idec];
            GoodChannel[ipro][idec] = true;
        }
    }
}

double KappaWidth(KAPPAS input)
{
    return pow(input.kc,2)*brc_SM + pow(input.kb,2)*brb_SM + pow(input.kmuon,2)*brmuon_SM + pow(input.ktau,2)*brtau_SM + pow(input.kw,2)*brW_SM + pow(input.kz,2)*brZ_SM + pow(input.kg,2)*brg_SM + pow(input.kga,2)*brgaga_SM + pow(input.kzga,2)*brZga_SM + brLeft_SM;
}

double VBF8Kappa(KAPPAS input)
{
// 1.210 from W-fusion
// 0.417 from Z-fusion
// The HEPfit data;
    return sqrt((pow(input.kw,2)*1.210 + pow(input.kz,2)*0.417)/(1.210+0.417));
}
double VBF13Kappa(KAPPAS input)
{
// 1.210 from W-fusion
// 0.417 from Z-fusion
// The HEPfit data; No 13 TeV directly data, using 8TeV instead
    return sqrt((pow(input.kw,2)*1.210 + pow(input.kz,2)*0.417)/(1.210+0.417));
}
double VH13Kappa(KAPPAS input)
{
// 0.8824 pb for ZH in SM
// 1.3690 pb for WH in SM
// Values are taken from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
    return sqrt((pow(input.kw,2)*1.369 + pow(input.kz,2)*0.8824)/(1.369+0.8824));
}

double GetKappa(KAPPAS input, int ID)
{
    switch(ID){
        case iU: return input.ku;
        case iD: return input.kd;
        case iC: return input.kc;
        case iS: return input.ks;
        case iT: return input.kt;
        case iB: return input.kb;
        case iE: return input.kele;
        case iM: return input.kmuon;
        case iL: return input.ktau;
        case iG: return input.kg;
        case iA: return input.kga;
        case iZ: return input.kz;
        case iW: return input.kw;
        case iZA: return input.kzga;
        case iVBF8: return VBF8Kappa(input);
        case iVBF13: return VBF13Kappa(input);
        case iVH13: return VH13Kappa(input);
        case iWid: return KappaWidth(input);
        default: return input.kz;
    }
}

STU_EXP::STU_EXP()
{
    SetUpExpData();
}
void STU_EXP::SetUpExpData()
{
    NObs = 0;
}
int STU_EXP::GetChiSquare(double S, double T, double U, double &chisquare)
{
    chisquare = 0.0;
    DOF = NObs;
    double Predictions[3] = {S,T,U};
    for (int i = 0; i < NObs; ++i)
    {
        for (int j = 0; j < NObs; ++j)
        {
            chisquare += (Predictions[i]-CentralValue[i])*SigmaInv[i][j]*(Predictions[j]-CentralValue[j]);
        }
    }
    return DOF;
}
STU_CEPC::STU_CEPC()
{
    SetUpExpData();
}
void STU_CEPC::SetUpExpData()
{
    NObs = 3;
    CentralValue[0] = 0.0;  Sigma[0] = 0.0246;
    CentralValue[1] = 0.0;  Sigma[1] = 0.0255;
    CentralValue[2] = 0.0;  Sigma[2] = 0.0208;

    Rho[0][0] = 1.000;  Rho[0][1] = 0.862;  Rho[0][2] =-0.373;
    Rho[1][0] = 0.862;  Rho[1][1] = 1.000;  Rho[1][2] =-0.735;
    Rho[2][0] =-0.373;  Rho[2][1] =-0.735;  Rho[2][2] = 1.000; 

    double det = Rho[0][2]*(Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1]) + Rho[0][1]*(Rho[1][0]*Rho[2][2]-Rho[1][2]*Rho[2][0]) + Rho[0][0]*(Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2]);

    SigmaInv[0][0] = (Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2])/Sigma[0]/Sigma[0]/det; SigmaInv[0][1] = (Rho[0][1]*Rho[2][2]-Rho[0][2]*Rho[2][1])/Sigma[0]/Sigma[1]/det; SigmaInv[0][2] = (Rho[0][2]*Rho[1][1]-Rho[0][1]*Rho[1][2])/Sigma[0]/Sigma[2]/det;
    SigmaInv[1][0] = (Rho[1][0]*Rho[2][2] - Rho[1][2]*Rho[2][0])/Sigma[0]/Sigma[1]/det; SigmaInv[1][1] = (Rho[0][2]*Rho[2][0]-Rho[0][0]*Rho[2][2])/Sigma[1]/Sigma[1]/det; SigmaInv[1][2] = (Rho[0][0]*Rho[1][2]-Rho[0][2]*Rho[1][0])/Sigma[1]/Sigma[2]/det;
    SigmaInv[2][0] = (Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1])/Sigma[0]/Sigma[2]/det; SigmaInv[2][1] = (Rho[0][0]*Rho[2][1]-Rho[0][1]*Rho[2][0])/Sigma[1]/Sigma[2]/det; SigmaInv[2][2] = (Rho[0][1]*Rho[1][0]-Rho[0][0]*Rho[1][1])/Sigma[2]/Sigma[2]/det;

}

STU_LHC::STU_LHC()
{
    SetUpExpData();
}
void STU_LHC::SetUpExpData()
{
    NObs = 3;
    CentralValue[0] = 0.04;  Sigma[0] = 0.11;
    CentralValue[1] = 0.09;  Sigma[1] = 0.14;
    CentralValue[2] = -0.02; Sigma[2] = 0.11;

    Rho[0][0] = 1.000;  Rho[0][1] = 0.920;  Rho[0][2] =-0.680;
    Rho[1][0] = 0.920;  Rho[1][1] = 1.000;  Rho[1][2] =-0.870;
    Rho[2][0] =-0.680;  Rho[2][1] =-0.870;  Rho[2][2] = 1.000; 

    double det = Rho[0][2]*(Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1]) + Rho[0][1]*(Rho[1][0]*Rho[2][2]-Rho[1][2]*Rho[2][0]) + Rho[0][0]*(Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2]);

    SigmaInv[0][0] = (Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2])/Sigma[0]/Sigma[0]/det; SigmaInv[0][1] = (Rho[0][1]*Rho[2][2]-Rho[0][2]*Rho[2][1])/Sigma[0]/Sigma[1]/det; SigmaInv[0][2] = (Rho[0][2]*Rho[1][1]-Rho[0][1]*Rho[1][2])/Sigma[0]/Sigma[2]/det;
    SigmaInv[1][0] = (Rho[1][0]*Rho[2][2] - Rho[1][2]*Rho[2][0])/Sigma[0]/Sigma[1]/det; SigmaInv[1][1] = (Rho[0][2]*Rho[2][0]-Rho[0][0]*Rho[2][2])/Sigma[1]/Sigma[1]/det; SigmaInv[1][2] = (Rho[0][0]*Rho[1][2]-Rho[0][2]*Rho[1][0])/Sigma[1]/Sigma[2]/det;
    SigmaInv[2][0] = (Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1])/Sigma[0]/Sigma[2]/det; SigmaInv[2][1] = (Rho[0][0]*Rho[2][1]-Rho[0][1]*Rho[2][0])/Sigma[1]/Sigma[2]/det; SigmaInv[2][2] = (Rho[0][1]*Rho[1][0]-Rho[0][0]*Rho[1][1])/Sigma[2]/Sigma[2]/det;

}

STU_ILC::STU_ILC()
{
    SetUpExpData();
}
void STU_ILC::SetUpExpData()
{
    NObs = 3;
    CentralValue[0] = 0.00;  Sigma[0] = 0.0352;
    CentralValue[1] = 0.00;  Sigma[1] = 0.0489;
    CentralValue[2] = 0.00;  Sigma[2] = 0.0376;

    Rho[0][0] = 1.000;  Rho[0][1] = 0.988;  Rho[0][2] =-0.879;
    Rho[1][0] = 0.988;  Rho[1][1] = 1.000;  Rho[1][2] =-0.909;
    Rho[2][0] =-0.879;  Rho[2][1] =-0.909;  Rho[2][2] = 1.000; 

    double det = Rho[0][2]*(Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1]) + Rho[0][1]*(Rho[1][0]*Rho[2][2]-Rho[1][2]*Rho[2][0]) + Rho[0][0]*(Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2]);

    SigmaInv[0][0] = (Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2])/Sigma[0]/Sigma[0]/det; SigmaInv[0][1] = (Rho[0][1]*Rho[2][2]-Rho[0][2]*Rho[2][1])/Sigma[0]/Sigma[1]/det; SigmaInv[0][2] = (Rho[0][2]*Rho[1][1]-Rho[0][1]*Rho[1][2])/Sigma[0]/Sigma[2]/det;
    SigmaInv[1][0] = (Rho[1][0]*Rho[2][2] - Rho[1][2]*Rho[2][0])/Sigma[0]/Sigma[1]/det; SigmaInv[1][1] = (Rho[0][2]*Rho[2][0]-Rho[0][0]*Rho[2][2])/Sigma[1]/Sigma[1]/det; SigmaInv[1][2] = (Rho[0][0]*Rho[1][2]-Rho[0][2]*Rho[1][0])/Sigma[1]/Sigma[2]/det;
    SigmaInv[2][0] = (Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1])/Sigma[0]/Sigma[2]/det; SigmaInv[2][1] = (Rho[0][0]*Rho[2][1]-Rho[0][1]*Rho[2][0])/Sigma[1]/Sigma[2]/det; SigmaInv[2][2] = (Rho[0][1]*Rho[1][0]-Rho[0][0]*Rho[1][1])/Sigma[2]/Sigma[2]/det;

}

STU_FCC::STU_FCC()
{
    SetUpExpData();
}
void STU_FCC::SetUpExpData()
{
    NObs = 3;
    CentralValue[0] = 0.00;  Sigma[0] = 0.0067;
    CentralValue[1] = 0.00;  Sigma[1] = 0.0053;
    CentralValue[2] = 0.00;  Sigma[2] = 0.0240;

    Rho[0][0] = 1.000;  Rho[0][1] = 0.812;  Rho[0][2] = 0.001;
    Rho[1][0] = 0.812;  Rho[1][1] = 1.000;  Rho[1][2] =-0.097;
    Rho[2][0] = 0.001;  Rho[2][1] =-0.097;  Rho[2][2] = 1.000; 

    double det = Rho[0][2]*(Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1]) + Rho[0][1]*(Rho[1][0]*Rho[2][2]-Rho[1][2]*Rho[2][0]) + Rho[0][0]*(Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2]);

    SigmaInv[0][0] = (Rho[1][2]*Rho[2][1] - Rho[1][1]*Rho[2][2])/Sigma[0]/Sigma[0]/det; SigmaInv[0][1] = (Rho[0][1]*Rho[2][2]-Rho[0][2]*Rho[2][1])/Sigma[0]/Sigma[1]/det; SigmaInv[0][2] = (Rho[0][2]*Rho[1][1]-Rho[0][1]*Rho[1][2])/Sigma[0]/Sigma[2]/det;
    SigmaInv[1][0] = (Rho[1][0]*Rho[2][2] - Rho[1][2]*Rho[2][0])/Sigma[0]/Sigma[1]/det; SigmaInv[1][1] = (Rho[0][2]*Rho[2][0]-Rho[0][0]*Rho[2][2])/Sigma[1]/Sigma[1]/det; SigmaInv[1][2] = (Rho[0][0]*Rho[1][2]-Rho[0][2]*Rho[1][0])/Sigma[1]/Sigma[2]/det;
    SigmaInv[2][0] = (Rho[1][1]*Rho[2][0] - Rho[1][0]*Rho[2][1])/Sigma[0]/Sigma[2]/det; SigmaInv[2][1] = (Rho[0][0]*Rho[2][1]-Rho[0][1]*Rho[2][0])/Sigma[1]/Sigma[2]/det; SigmaInv[2][2] = (Rho[0][1]*Rho[1][0]-Rho[0][0]*Rho[1][1])/Sigma[2]/Sigma[2]/det;

}

mu_LHC8TeV::mu_LHC8TeV()
{
    SetUpExpData();
}
void mu_LHC8TeV::SetUpExpData()
{
// From 1606.02266
    NProd = 5;
    NDecay = 5;

    ProdID[0] = iG;
    ProdID[1] = iVBF8;
    ProdID[2] = iW;
    ProdID[3] = iZ;
    ProdID[4] = iT;

    DecayID[0] = iA;
    DecayID[1] = iW;
    DecayID[2] = iZ;
    DecayID[3] = iL;
    DecayID[4] = iB;

    CentralValue[0][0] = 1.10; GoodChannel[0][0] = true;  UPError[0][0] = 0.23; DOError[0][0] = 0.22; //ggF > aa
    CentralValue[1][0] = 1.30; GoodChannel[1][0] = true;  UPError[1][0] = 0.50; DOError[1][0] = 0.50; //VBF > aa
    CentralValue[2][0] = 0.50; GoodChannel[2][0] = true;  UPError[2][0] = 1.30; DOError[2][0] = 1.20; //Wh > aa
    CentralValue[3][0] = 0.50; GoodChannel[3][0] = true;  UPError[3][0] = 3.00; DOError[3][0] = 2.50; //Zh > aa
    CentralValue[4][0] = 2.20; GoodChannel[4][0] = true;  UPError[4][0] = 1.60; DOError[4][0] = 1.30; //tth > aa

    CentralValue[0][1] = 0.84; GoodChannel[0][1] = true;  UPError[0][1] = 0.17; DOError[0][1] = 0.17; //ggF > WW
    CentralValue[1][1] = 1.20; GoodChannel[1][1] = true;  UPError[1][1] = 0.40; DOError[1][1] = 0.40; //VBF > WW
    CentralValue[2][1] = 1.60; GoodChannel[2][1] = true;  UPError[2][1] = 1.20; DOError[2][1] = 1.00; //Wh > WW
    CentralValue[3][1] = 5.90; GoodChannel[3][1] = true;  UPError[3][1] = 2.60; DOError[3][1] = 2.20; //Zh > WW
    CentralValue[4][1] = 5.00; GoodChannel[4][1] = true;  UPError[4][1] = 1.80; DOError[4][1] = 1.70; //tth > WW

    CentralValue[0][2] = 1.13; GoodChannel[0][2] = true;  UPError[0][2] = 0.34; DOError[0][2] = 0.31; //ggF > ZZ
    CentralValue[1][2] = 0.10; GoodChannel[1][2] = true;  UPError[1][2] = 1.10; DOError[1][2] = 0.60; //VBF > ZZ
    CentralValue[2][2] = 1.00; GoodChannel[2][2] = false; UPError[2][2] = 1.00; DOError[2][2] = 1.00; //Wh > ZZ
    CentralValue[3][2] = 1.00; GoodChannel[3][2] = false; UPError[3][2] = 1.00; DOError[3][2] = 1.00; //Zh > ZZ
    CentralValue[4][2] = 1.00; GoodChannel[4][2] = false; UPError[4][2] = 1.00; DOError[4][2] = 1.00; //tth > ZZ

    CentralValue[0][3] = 1.00; GoodChannel[0][3] = true;  UPError[0][3] = 0.60; DOError[0][3] = 0.60; //ggF > ta ta
    CentralValue[1][3] = 1.30; GoodChannel[1][3] = true;  UPError[1][3] = 0.40; DOError[1][3] = 0.40; //VBF > ta ta
    CentralValue[2][3] = -1.4; GoodChannel[2][3] = true;  UPError[2][3] = 1.40; DOError[2][3] = 1.40; //Wh > ta ta
    CentralValue[3][3] = 2.20; GoodChannel[3][3] = true;  UPError[3][3] = 2.20; DOError[3][3] = 1.80; //Zh > ta ta
    CentralValue[4][3] = -1.9; GoodChannel[4][3] = true;  UPError[4][3] = 3.70; DOError[4][3] = 3.30; //tth > ta ta

    CentralValue[0][4] = 1.00; GoodChannel[0][4] = false; UPError[0][4] = 1.00; DOError[0][4] = 1.00; //ggF > bb
    CentralValue[1][4] = 1.00; GoodChannel[1][4] = false; UPError[1][4] = 1.00; DOError[1][4] = 1.00; //VBF > bb
    CentralValue[2][4] = 1.00; GoodChannel[2][4] = true;  UPError[2][4] = 0.50; DOError[2][4] = 0.50; //Wh > bb
    CentralValue[3][4] = 0.40; GoodChannel[3][4] = true;  UPError[3][4] = 0.40; DOError[3][4] = 0.40; //Zh > bb
    CentralValue[4][4] = 1.10; GoodChannel[4][4] = true;  UPError[4][4] = 1.00; DOError[4][4] = 1.00; //tth > bb
}

mu_ATLAS13TeV::mu_ATLAS13TeV()
{
    SetUpExpData();
}
void mu_ATLAS13TeV::SetUpExpData()
{
// From ATLAS-CONF-2018-031 and ATLAS-CONF-2016-112
    NProd = 4;
    NDecay = 5;

    ProdID[0] = iG;
    ProdID[1] = iVBF13;
    ProdID[2] = iVH13;
    ProdID[3] = iT;

    DecayID[0] = iA;
    DecayID[1] = iZ;
    DecayID[2] = iW;
    DecayID[3] = iB;
    DecayID[4] = iL;

    CentralValue[0][0] = 0.97; GoodChannel[0][0] = true;  UPError[0][0] = 0.14; DOError[0][0] = 0.14; //ggF > aa
    CentralValue[1][0] = 1.42; GoodChannel[1][0] = true;  UPError[1][0] = 0.43; DOError[1][0] = 0.37; //VBF > aa
    CentralValue[2][0] = 1.09; GoodChannel[2][0] = true;  UPError[2][0] = 0.59; DOError[2][0] = 0.55; //Vh > aa
    CentralValue[3][0] = 1.14; GoodChannel[3][0] = true;  UPError[3][0] = 0.43; DOError[3][0] = 0.37; //tth > aa

    CentralValue[0][1] = 1.04; GoodChannel[0][1] = true;  UPError[0][1] = 0.17; DOError[0][1] = 0.15; //ggF > ZZ
    CentralValue[1][1] = 2.96; GoodChannel[1][1] = true;  UPError[1][1] = 0.97; DOError[1][1] = 0.83; //VBF > ZZ
    CentralValue[2][1] = 0.70; GoodChannel[2][1] = true;  UPError[2][1] = 1.22; DOError[2][1] = 0.79; //Vh > ZZ
    CentralValue[3][1] = 1.49; GoodChannel[3][1] = true;  UPError[3][1] = 0.60; DOError[3][1] = 0.57; //tth > ZZ

    CentralValue[0][2] = 1.20; GoodChannel[0][2] = true;  UPError[0][2] = 0.20; DOError[0][2] = 0.19; //ggF > WW
    CentralValue[1][2] = 0.61; GoodChannel[1][2] = true;  UPError[1][2] = 0.37; DOError[1][2] = 0.36; //VBF > WW
    CentralValue[2][2] = 3.20; GoodChannel[2][2] = true;  UPError[2][2] = 4.40; DOError[2][2] = 4.20; //Vh > WW
    CentralValue[3][2] = 1.00; GoodChannel[3][2] = false; UPError[3][2] = 1.00; DOError[3][2] = 1.00; //tth > WW

    CentralValue[0][3] = 1.00; GoodChannel[0][3] = false; UPError[0][3] = 1.00; DOError[0][3] = 1.00; //ggF > bb
    CentralValue[1][3] = 1.00; GoodChannel[1][3] = false; UPError[1][3] = 1.00; DOError[1][3] = 1.00; //VBF > bb
    CentralValue[2][3] = 1.19; GoodChannel[2][3] = true;  UPError[2][3] = 0.41; DOError[2][3] = 0.36; //Vh > bb
    CentralValue[3][3] = 0.80; GoodChannel[3][3] = true;  UPError[3][3] = 0.59; DOError[3][3] = 0.60; //tth > bb

    CentralValue[0][4] = 0.97; GoodChannel[0][4] = true;  UPError[0][4] = 0.59; DOError[0][4] = 0.51; //ggF > tata
    CentralValue[1][4] = 1.16; GoodChannel[1][4] = true;  UPError[1][4] = 0.57; DOError[1][4] = 0.42; //VBF > tata
    CentralValue[2][4] = 1.00; GoodChannel[2][4] = false; UPError[2][4] = 1.00; DOError[2][4] = 1.00; //Vh > tata
    CentralValue[3][4] = 1.41; GoodChannel[3][4] = true;  UPError[3][4] = 1.13; DOError[3][4] = 0.97; //tth > tata
}

mu_CMS13TeV::mu_CMS13TeV()
{
    SetUpExpData();
}
void mu_CMS13TeV::SetUpExpData()
{
// From 1809.10733 Table. 3
    NProd = 5;
    NDecay = 6;

    ProdID[0] = iG;
    ProdID[1] = iVBF13;
    ProdID[2] = iW;
    ProdID[3] = iZ;
    ProdID[4] = iT;

    DecayID[0] = iA;
    DecayID[1] = iZ;
    DecayID[2] = iW;
    DecayID[3] = iL;
    DecayID[4] = iB;
    DecayID[5] = iM;

    CentralValue[0][0] = 1.16; GoodChannel[0][0] = true;  UPError[0][0] = 0.21; DOError[0][0] = 0.18; //ggF > aa
    CentralValue[1][0] = 0.67; GoodChannel[1][0] = true;  UPError[1][0] = 0.59; DOError[1][0] = 0.46; //VBF > aa
    CentralValue[2][0] = 3.76; GoodChannel[2][0] = true;  UPError[2][0] = 1.48; DOError[2][0] = 1.35; //Wh > aa
    CentralValue[3][0] = 0.00; GoodChannel[3][0] = true;  UPError[3][0] = 1.14; DOError[3][0] = 0.01; //Zh > aa
    CentralValue[4][0] = 2.18; GoodChannel[4][0] = true;  UPError[4][0] = 0.88; DOError[4][0] = 0.75; //tth > aa

    CentralValue[0][1] = 1.22; GoodChannel[0][1] = true;  UPError[0][1] = 0.23; DOError[0][1] = 0.21; //ggF > ZZ
    CentralValue[1][1] =-0.09; GoodChannel[1][1] = true;  UPError[1][1] = 1.02; DOError[1][1] = 0.76; //VBF > ZZ
    CentralValue[2][1] = 0.00; GoodChannel[2][1] = true;  UPError[2][1] = 2.33; DOError[2][1] = 0.01; //Wh > ZZ
    CentralValue[3][1] = 0.00; GoodChannel[3][1] = true;  UPError[3][1] = 4.26; DOError[3][1] = 0.01; //Zh > ZZ
    CentralValue[4][1] = 0.00; GoodChannel[4][1] = true;  UPError[4][1] = 1.50; DOError[4][1] = 0.01; //tth > ZZ

    CentralValue[0][2] = 1.35; GoodChannel[0][2] = true;  UPError[0][2] = 0.21; DOError[0][2] = 0.19; //ggF > WW
    CentralValue[1][2] = 0.28; GoodChannel[1][2] = true;  UPError[1][2] = 0.64; DOError[1][2] = 0.60; //VBF > WW
    CentralValue[2][2] = 3.91; GoodChannel[2][2] = true;  UPError[2][2] = 2.26; DOError[2][2] = 2.01; //Wh > WW
    CentralValue[3][2] = 0.96; GoodChannel[3][2] = true;  UPError[3][2] = 1.81; DOError[3][2] = 1.46; //Zh > WW
    CentralValue[4][2] = 1.60; GoodChannel[4][2] = true;  UPError[4][2] = 0.65; DOError[4][2] = 0.59; //tth > WW

    CentralValue[0][3] = 1.05; GoodChannel[0][3] = true;  UPError[0][3] = 0.53; DOError[0][3] = 0.47; //ggF > ta ta
    CentralValue[1][3] = 1.12; GoodChannel[1][3] = true;  UPError[1][3] = 0.45; DOError[1][3] = 0.43; //VBF > ta ta
    CentralValue[2][3] = 1.00; GoodChannel[2][3] = false; UPError[2][3] = 1.00; DOError[2][3] = 1.00; //Wh > ta ta
    CentralValue[3][3] = 1.00; GoodChannel[3][3] = false; UPError[3][3] = 1.00; DOError[3][3] = 1.00; //Zh > ta ta
    CentralValue[4][3] = 0.23; GoodChannel[4][3] = true;  UPError[4][3] = 1.03; DOError[4][3] = 0.88; //tth > ta ta

    CentralValue[0][4] = 2.51; GoodChannel[0][4] = true;  UPError[0][4] = 2.43; DOError[0][4] = 2.01; //ggF > bb
    CentralValue[1][4] = 1.00; GoodChannel[1][4] = false; UPError[1][4] = 1.00; DOError[1][4] = 1.00; //VBF > bb
    CentralValue[2][4] = 1.73; GoodChannel[2][4] = true;  UPError[2][4] = 0.70; DOError[2][4] = 0.68; //Wh > bb
    CentralValue[3][4] = 0.99; GoodChannel[3][4] = true;  UPError[3][4] = 0.47; DOError[3][4] = 0.45; //Zh > bb
    CentralValue[4][4] = 0.91; GoodChannel[4][4] = true;  UPError[4][4] = 0.45; DOError[4][4] = 0.43; //tth > bb

    CentralValue[0][5] = 0.31; GoodChannel[0][5] = true;  UPError[0][5] = 1.80; DOError[0][5] = 1.79; //ggF > mumu
    CentralValue[1][5] = 2.72; GoodChannel[1][5] = true;  UPError[1][5] = 7.12; DOError[1][5] = 7.03; //VBF > mumu
    CentralValue[2][5] = 1.00; GoodChannel[2][5] = false; UPError[2][5] = 1.00; DOError[2][5] = 1.00; //Wh > mumu
    CentralValue[3][5] = 1.00; GoodChannel[3][5] = false; UPError[3][5] = 1.00; DOError[3][5] = 1.00; //Zh > mumu
    CentralValue[4][5] = 1.00; GoodChannel[4][5] = false; UPError[4][5] = 1.00; DOError[4][5] = 1.00; //tth > mumu

}

mu_HLLHC300::mu_HLLHC300()
{
    SetUpExpData();
}
void mu_HLLHC300::SetUpExpData()
{
    NDecay = 7; // gamma-gamma ZZ WW Z-gamma bb tata mumu
    NProd = 5; // ggF VBF Wh Zh tth
    ProdID[0] = iG;
    ProdID[1] = iVBF13;
    ProdID[2] = iW;
    ProdID[3] = iZ;
    ProdID[4] = iT;

    DecayID[0] = iA;
    DecayID[1] = iZ;
    DecayID[2] = iW;
    DecayID[3] = iZA;
    DecayID[4] = iB;
    DecayID[5] = iL;
    DecayID[6] = iM;

    NChannel[0] = 7; // For gamma-gamma, 7 channels
    NChannel[1] = 5; // For ZZ, 5 channels
    NChannel[2] = 4; // For WW, 4 channels
    NChannel[3] = 1; // For Z-gamma, 1 channel
    NChannel[4] = 3; // For bb, 3 channels
    NChannel[5] = 1; // For tata, 1 channel
    NChannel[6] = 3; // For mumu, 3 channels


// First index for channel, second index for decay 
    UPError[0][0] = 0.09; // comb  gaga
    UPError[1][0] = 0.12; // 0j gaga
    UPError[2][0] = 0.14; // 1j gaga
    UPError[3][0] = 0.43; // VBF-like gaga
    UPError[4][0] = 0.48; // WH-like gaga
    UPError[5][0] = 0.85; // ZH-like gaga
    UPError[6][0] = 0.36; // ttH-like gaga

    UPError[0][1] = 0.07; // comb ZZ
    UPError[1][1] = 0.34; // VH-like ZZ
    UPError[2][1] = 0.48; // ttH-like ZZ
    UPError[3][1] = 0.33; // VBF-like ZZ
    UPError[4][1] = 0.07; // ggF-like ZZ

    UPError[0][2] = 0.08; // comb WW
    UPError[1][2] = 0.09; // 0j WW
    UPError[2][2] = 0.18; // 1j WW
    UPError[3][2] = 0.20; // VBF-like WW

    UPError[0][3] = 0.44; // incl. Z-ga

    UPError[0][4] = 0.26; // comb bb
    UPError[1][4] = 0.56; // WH-like bb
    UPError[2][4] = 0.29; // ZH-like bb

    UPError[0][5] = 0.18; // VBF-like tata

    UPError[0][6] = 0.38; // comb mumu
    UPError[1][6] = 0.45; // incl. mumu
    UPError[2][6] = 0.72; // ttH-like mumu

    for (int idec = 0; idec < NDecay; ++idec)
    {
        for (int ichan = 0; ichan < NChannel[idec]; ++ichan)
        {
            DOError[ichan][idec] = UPError[ichan][idec];
            CentralValue[ichan][idec] = 1.0;
            for (int ipro = 0; ipro < NProd; ++ipro)
            {
                Contamination[ipro][ichan][idec] = 0.0;
            }
        }
    }

// First channel for Production mode, second for channel, third for decay
    Contamination[0][0][0] = 49.85/56.92;
    Contamination[1][0][0] = 4.18/56.92;
    Contamination[2][0][0] = 1.5/56.92;
    Contamination[3][0][0] = 0.88/56.92;
    Contamination[4][0][0] = 0.61/56.92;

    Contamination[0][1][0] = 97.0/100.0;
    Contamination[1][1][0] = 3.0/100.0;

    Contamination[0][2][0] = 86.0/100.0;
    Contamination[1][2][0] = 14.0/100.0;

    Contamination[0][3][0] = 0.3;
    Contamination[1][3][0] = 0.7;

    Contamination[2][4][0] = 72.0/77.0;
    Contamination[4][4][0] = 5.0/77.0;

    Contamination[2][5][0] = 10.0/102.0;
    Contamination[3][5][0] = 88.0/102.0;
    Contamination[4][5][0] = 4.0/102.0;

    Contamination[2][6][0] = 7.0/206.0;
    Contamination[3][6][0] = 12.0/206.0;
    Contamination[4][6][0] = 187.0/206.0;

    Contamination[0][0][1] = 49.85/56.92;
    Contamination[1][0][1] = 4.18/56.92;
    Contamination[2][0][1] = 1.5/56.92;
    Contamination[3][0][1] = 0.88/56.92;
    Contamination[4][0][1] = 0.61/56.92;

    Contamination[0][1][1] = 22.0/72.5;
    Contamination[1][1][1] = 6.6/72.5;
    Contamination[2][1][1] = 25.0/72.5;
    Contamination[3][1][1] = 8.8/72.5;
    Contamination[4][1][1] = 10.1/72.5;

    Contamination[0][2][1] = 3.1/35.4;
    Contamination[1][2][1] = 0.6/35.4;
    Contamination[2][2][1] = 0.6/35.4;
    Contamination[3][2][1] = 1.1/35.4;
    Contamination[4][2][1] = 30.0/35.4;

    Contamination[0][3][1] = 41.0/97.1;
    Contamination[1][3][1] = 54.0/97.1;
    Contamination[2][3][1] = 0.7/97.1;
    Contamination[3][3][1] = 0.4/97.1;
    Contamination[4][3][1] = 1.0/97.1;

    Contamination[0][4][1] = 3380.0/3809.0;
    Contamination[1][4][1] = 274.0/3809.0;
    Contamination[2][4][1] = 77.0/3809.0;
    Contamination[3][4][1] = 53.0/3809.0;
    Contamination[4][4][1] = 25.0/3809.0;

    Contamination[0][0][2] = 49.85/56.92;
    Contamination[1][0][2] = 4.18/56.92;
    Contamination[2][0][2] = 1.5/56.92;
    Contamination[3][0][2] = 0.88/56.92;
    Contamination[4][0][2] = 0.61/56.92;

    Contamination[0][1][2] = 4085.0/4184.0;
    Contamination[1][1][2] = 99.0/4184.0;

    Contamination[0][2][2] = 20050.0/22375.0;
    Contamination[1][2][2] = 2325.0/22375.0;

    Contamination[0][3][2] = 9.0/59.0;
    Contamination[1][3][2] = 50.0/59.0;

    Contamination[0][0][3] = 12.79/28.13;
    Contamination[1][0][3] = 15.34/28.13;

    Contamination[0][0][4] = 49.85/56.92;
    Contamination[1][0][4] = 4.18/56.92;
    Contamination[2][0][4] = 1.5/56.92;
    Contamination[3][0][4] = 0.88/56.92;
    Contamination[4][0][4] = 0.61/56.92;

    Contamination[0][1][4] = 1.0;

    Contamination[3][2][4] = 560.0/636.0;
    Contamination[4][2][4] = 76.0/636.0;

    Contamination[0][0][5] = 1641.0/2538.0;
    Contamination[1][0][5] = 897.0/2538.0;

    Contamination[0][0][6] = 49.85/56.92;
    Contamination[1][0][6] = 4.18/56.92;
    Contamination[2][0][6] = 1.5/56.92;
    Contamination[3][0][6] = 0.88/56.92;
    Contamination[4][0][6] = 0.61/56.92;

    Contamination[0][1][6] = 1510.0/1725.0;
    Contamination[1][1][6] = 125.0/1725.0;
    Contamination[2][1][6] = 45.0/1725.0;
    Contamination[3][1][6] = 27.0/1725.0;
    Contamination[4][1][6] = 18.0/1725.0;

    Contamination[4][2][6] = 1.0;

}
int mu_HLLHC300::GetChiSquare(KAPPAS input, double &chisquare)
{
    chisquare = 0.0;
    DOF = 0;
    double kwid = GetKappa(input,iWid);
    double kprod2;
    double kdecay;
    double mu;
    for (int idec = 0; idec < NDecay; ++idec)
    {
        for (int ichan = 0; ichan < NChannel[idec]; ++ichan)
        {
            kprod2 = 0.0;
            for (int ipro = 0; ipro < NProd; ++ipro)
            {
                kprod2 += Contamination[ipro][ichan][idec]*pow(GetKappa(input,ProdID[ipro]),2);
            }
            kdecay = GetKappa(input,DecayID[idec]);
            mu = kprod2*pow(kdecay,2)/kwid;
            DOF+=1;
            chisquare += pow((mu-CentralValue[ichan][idec])/UPError[ichan][idec],2);
        }
    }
    return DOF;
}

mu_HLLHC3000::mu_HLLHC3000()
{
    SetUpExpData();
}
void mu_HLLHC3000::SetUpExpData()
{
    NDecay = 7; // gamma-gamma ZZ WW Z-gamma bb tata mumu
    NProd = 5; // ggF VBF Wh Zh tth
    ProdID[0] = iG;
    ProdID[1] = iVBF13;
    ProdID[2] = iW;
    ProdID[3] = iZ;
    ProdID[4] = iT;

    DecayID[0] = iA;
    DecayID[1] = iZ;
    DecayID[2] = iW;
    DecayID[3] = iZA;
    DecayID[4] = iB;
    DecayID[5] = iL;
    DecayID[6] = iM;

    NChannel[0] = 7; // For gamma-gamma, 7 channels
    NChannel[1] = 5; // For ZZ, 5 channels
    NChannel[2] = 4; // For WW, 4 channels
    NChannel[3] = 1; // For Z-gamma, 1 channel
    NChannel[4] = 3; // For bb, 3 channels
    NChannel[5] = 1; // For tata, 1 channel
    NChannel[6] = 3; // For mumu, 3 channels


// First index for channel, second index for decay 
    UPError[0][0] = 0.04; // comb  gaga
    UPError[1][0] = 0.05; // 0j gaga
    UPError[2][0] = 0.23; // 1j gaga
    UPError[3][0] = 0.15; // VBF-like gaga
    UPError[4][0] = 0.17; // WH-like gaga
    UPError[5][0] = 0.28; // ZH-like gaga
    UPError[6][0] = 0.12; // ttH-like gaga

    UPError[0][1] = 0.04; // comb ZZ
    UPError[1][1] = 0.12; // VH-like ZZ
    UPError[2][1] = 0.20; // ttH-like ZZ
    UPError[3][1] = 0.16; // VBF-like ZZ
    UPError[4][1] = 0.04; // ggF-like ZZ

    UPError[0][2] = 0.05; // comb WW
    UPError[1][2] = 0.05; // 0j WW
    UPError[2][2] = 0.10; // 1j WW
    UPError[3][2] = 0.09; // VBF-like WW

    UPError[0][3] = 0.27; // incl. Z-ga

    UPError[0][4] = 0.12; // comb bb
    UPError[1][4] = 0.36; // WH-like bb
    UPError[2][4] = 0.13; // ZH-like bb

    UPError[0][5] = 0.15; // VBF-like tata

    UPError[0][6] = 0.12; // comb mumu
    UPError[1][6] = 0.14; // incl. mumu
    UPError[2][6] = 0.23; // ttH-like mumu

    for (int idec = 0; idec < NDecay; ++idec)
    {
        for (int ichan = 0; ichan < NChannel[idec]; ++ichan)
        {
            DOError[ichan][idec] = UPError[ichan][idec];
            CentralValue[ichan][idec] = 1.0;
            for (int ipro = 0; ipro < NProd; ++ipro)
            {
                Contamination[ipro][ichan][idec] = 0.0;
            }
        }
    }

// First channel for Production mode, second for channel, third for decay
    Contamination[0][0][0] = 49.85/56.92;
    Contamination[1][0][0] = 4.18/56.92;
    Contamination[2][0][0] = 1.5/56.92;
    Contamination[3][0][0] = 0.88/56.92;
    Contamination[4][0][0] = 0.61/56.92;

    Contamination[0][1][0] = 97.0/100.0;
    Contamination[1][1][0] = 3.0/100.0;

    Contamination[0][2][0] = 86.0/100.0;
    Contamination[1][2][0] = 14.0/100.0;

    Contamination[0][3][0] = 0.4;
    Contamination[1][3][0] = 0.6;

    Contamination[2][4][0] = 72.0/77.0;
    Contamination[4][4][0] = 5.0/77.0;

    Contamination[2][5][0] = 10.0/102.0;
    Contamination[3][5][0] = 88.0/102.0;
    Contamination[4][5][0] = 4.0/102.0;

    Contamination[2][6][0] = 7.0/206.0;
    Contamination[3][6][0] = 12.0/206.0;
    Contamination[4][6][0] = 187.0/206.0;

    Contamination[0][0][1] = 49.85/56.92;
    Contamination[1][0][1] = 4.18/56.92;
    Contamination[2][0][1] = 1.5/56.92;
    Contamination[3][0][1] = 0.88/56.92;
    Contamination[4][0][1] = 0.61/56.92;

    Contamination[0][1][1] = 22.0/72.5;
    Contamination[1][1][1] = 6.6/72.5;
    Contamination[2][1][1] = 25.0/72.5;
    Contamination[3][1][1] = 8.8/72.5;
    Contamination[4][1][1] = 10.1/72.5;

    Contamination[0][2][1] = 3.1/35.4;
    Contamination[1][2][1] = 0.6/35.4;
    Contamination[2][2][1] = 0.6/35.4;
    Contamination[3][2][1] = 1.1/35.4;
    Contamination[4][2][1] = 30.0/35.4;

    Contamination[0][3][1] = 41.0/97.1;
    Contamination[1][3][1] = 54.0/97.1;
    Contamination[2][3][1] = 0.7/97.1;
    Contamination[3][3][1] = 0.4/97.1;
    Contamination[4][3][1] = 1.0/97.1;

    Contamination[0][4][1] = 3380.0/3809.0;
    Contamination[1][4][1] = 274.0/3809.0;
    Contamination[2][4][1] = 77.0/3809.0;
    Contamination[3][4][1] = 53.0/3809.0;
    Contamination[4][4][1] = 25.0/3809.0;

    Contamination[0][0][2] = 49.85/56.92;
    Contamination[1][0][2] = 4.18/56.92;
    Contamination[2][0][2] = 1.5/56.92;
    Contamination[3][0][2] = 0.88/56.92;
    Contamination[4][0][2] = 0.61/56.92;

    Contamination[0][1][2] = 4085.0/4184.0;
    Contamination[1][1][2] = 99.0/4184.0;

    Contamination[0][2][2] = 20050.0/22375.0;
    Contamination[1][2][2] = 2325.0/22375.0;

    Contamination[0][3][2] = 9.0/59.0;
    Contamination[1][3][2] = 50.0/59.0;

    Contamination[0][0][3] = 12.79/28.13;
    Contamination[1][0][3] = 15.34/28.13;

    Contamination[0][0][4] = 49.85/56.92;
    Contamination[1][0][4] = 4.18/56.92;
    Contamination[2][0][4] = 1.5/56.92;
    Contamination[3][0][4] = 0.88/56.92;
    Contamination[4][0][4] = 0.61/56.92;

    Contamination[0][1][4] = 1.0;

    Contamination[3][2][4] = 560.0/636.0;
    Contamination[4][2][4] = 76.0/636.0;

    Contamination[0][0][5] = 1641.0/2538.0;
    Contamination[1][0][5] = 897.0/2538.0;

    Contamination[0][0][6] = 49.85/56.92;
    Contamination[1][0][6] = 4.18/56.92;
    Contamination[2][0][6] = 1.5/56.92;
    Contamination[3][0][6] = 0.88/56.92;
    Contamination[4][0][6] = 0.61/56.92;

    Contamination[0][1][6] = 1510.0/1725.0;
    Contamination[1][1][6] = 125.0/1725.0;
    Contamination[2][1][6] = 45.0/1725.0;
    Contamination[3][1][6] = 27.0/1725.0;
    Contamination[4][1][6] = 18.0/1725.0;

    Contamination[4][2][6] = 1.0;

}
int mu_HLLHC3000::GetChiSquare(KAPPAS input, double &chisquare)
{
    chisquare = 0.0;
    DOF = 0;
    double kwid = GetKappa(input,iWid);
    double kprod2;
    double kdecay;
    double mu;
    for (int idec = 0; idec < NDecay; ++idec)
    {
        for (int ichan = 0; ichan < NChannel[idec]; ++ichan)
        {
            kprod2 = 0.0;
            for (int ipro = 0; ipro < NProd; ++ipro)
            {
                kprod2 += Contamination[ipro][ichan][idec]*pow(GetKappa(input,ProdID[ipro]),2);
            }
            kdecay = GetKappa(input,DecayID[idec]);
            mu = kprod2*pow(kdecay,2)/kwid;
            DOF+=1;
            chisquare += pow((mu-CentralValue[ichan][idec])/UPError[ichan][idec],2);
        }
    }
    return DOF;
}

double ChiSquare95Table[61] = {0,3.841,5.991,7.815,9.488,11.07,12.592,14.067,15.507,16.919,18.307,19.675,21.026,22.362,23.685,24.996,26.296,27.587,28.869,30.144,31.410,32.671,33.924,35.172,36.415,37.652,38.885,40.113,42.557,43.773,44.985,46.194,47.400,48.602,49.802,50.998,52.192,53.384,54.572,55.758,56.942,58.124,59.304,60.481,61.656,62.830,64.001,65.171,66.339,67.505,68.669,69.832,70.993,72.153,73.311,74.468,75.624,76.778,77.931,79.082};

void ChiSquare_Test(double chisquare, int DOF, bool &passed)
{
    passed = chisquare < ChiSquare95Table[DOF];
}
void HiggsSignalStrength_Test(int Exps, KAPPAS input, double &chi2mu, int &DOF, bool &passed)
{
    DOF = 0;
    chi2mu = 0.0;
    double chisqtemp;
    for (int i = 0; i < NEXPmu; ++i)
    {
        if ((Exps>>i) & 1)
        {
            DOF += AllmuExps[i].GetChiSquare(input,chisqtemp);
            chi2mu += chisqtemp;
        }
    }
    ChiSquare_Test(chi2mu,DOF,passed);
}
void STU_Test(int Exps, double S, double T, double U, double &chi2STU, int &DOF, bool &passed)
{
   DOF = 0;
   chi2STU = 0.0;
   double chisqtemp;
   for (int i = 0; i < NEXPSTU; ++i)
   {
       if ((Exps>>i) & 1)
       {
           DOF += AllSTUExps[i].GetChiSquare(S,T,U,chisqtemp);
           chi2STU += chisqtemp;
       }
   }
   ChiSquare_Test(chi2STU,DOF,passed);
}
void HiggsSignalStrengthSTU_Test(int muExps, int STUExps, KAPPAS input, double S, double T, double U, double &chi2mu, double &chi2STU, int &DOF, bool &passed)
{
    // DOF = 0;
    // // chisquare = 0.0;
    // // double chisqtemp;
    // DOF+=muExpCEPC.GetChiSquare(input,chi2mu);
    // // chisquare += chisqtemp;
    // DOF+=STUExpCEPC.GetChiSquare(S,T,U,chi2STU);
    // // chisquare += chisqtemp;
    // passed = (chi2mu+chi2STU) < ChiSquare95Table[DOF];
    int DOFSTU;
    int DOFmu;
    bool mupassed;
    bool STUpassed;
    HiggsSignalStrength_Test(muExps,input,chi2mu,DOFmu,mupassed);
    STU_Test(STUExps,S,T,U,chi2STU,DOFSTU,STUpassed);
    DOF = DOFSTU + DOFmu;
    ChiSquare_Test(chi2STU+chi2mu,DOF,passed);
}
