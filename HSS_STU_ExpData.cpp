#include "HSS_STU_ExpData.h"

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
    int DOF = 0;
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
        }
    }
    return DOF;
}

CEPC::CEPC()
{
    SetUpExpData();
}

void CEPC::SetUpExpData()
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
CEPC_STU::CEPC_STU()
{
    SetUpExpData();
}
void CEPC_STU::SetUpExpData()
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

double ChiSquare95Table[31] = {0,3.841,5.991,7.815,9.488,11.07,12.592,14.067,15.507,16.919,18.307,19.675,21.026,22.362,23.685,24.996,26.296,27.587,28.869,30.144,31.410,32.671,33.924,35.172,36.415,37.652,38.885,40.113,42.557,43.773};

void HiggsSignalStrengthSTU_Test(KAPPAS input, double S, double T, double U, double &chi2mu, double &chi2STU, int &DOF, bool &passed)
{
    DOF = 0;
    // chisquare = 0.0;
    // double chisqtemp;
    DOF+=mu_CEPC.GetChiSquare(input,chi2mu);
    // chisquare += chisqtemp;
    DOF+=STU_CEPC.GetChiSquare(S,T,U,chi2STU);
    // chisquare += chisqtemp;
    passed = (chi2mu+chi2STU) < ChiSquare95Table[DOF];
}
