///______________________________________________________________________________
///______________________________________________________________________________
///______________________________________________________________________________
///             
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestHitsKalmanFilterProcess.cxx
///
///             MAY 2019 : lujian
//
///_______________________________________________________________________________
#include <Math/SMatrix.h>
#include <Math/SVector.h>
#include <TLegend.h>
#include <TRandom3.h>
#include "TRestHitsKalmanFilterProcess.h"
#include <fstream>
#include <sstream>
using namespace std;

typedef ROOT::Math::SMatrix<double,6,6>  SMatrix66;
typedef ROOT::Math::SMatrix<double,3,6>  SMatrix36;
typedef ROOT::Math::SMatrix<double,6,3>  SMatrix63;
typedef ROOT::Math::SMatrix<double,3,3>  SMatrix33;
typedef ROOT::Math::SVector<double,6>  SVector6;
typedef ROOT::Math::SVector<double,3>  SVector3;
typedef ROOT::Math::SVector<double,3>  SVector3;
Double_t lightSpeed=300000000; // units m/s
//Double_t dZ=3.0; //units mm
Double_t MeasureError=3.0; //units mm

ClassImp(TRestHitsKalmanFilterProcess)
//______________________________________________________________________________
TRestHitsKalmanFilterProcess::TRestHitsKalmanFilterProcess( )
{
    Initialize();
}

//______________________________________________________________________________
TRestHitsKalmanFilterProcess::TRestHitsKalmanFilterProcess( char *cfgFileName )
{
    Initialize();

    if( LoadConfigFromFile( cfgFileName ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
TRestHitsKalmanFilterProcess::~TRestHitsKalmanFilterProcess( )
{
    delete fInputHitsEvent;
    delete fOutputHitsEvent;
}

void TRestHitsKalmanFilterProcess::LoadDefaultConfig( )
{
    SetName( "HitsKalmanFilter" );
    SetTitle( "Default config" );

  /*  fStartingDistance = 0.5;
    fMinimumDistance  = 3;
    fDistanceFactor   = 1.5;
    fMaxNodes         = 30;
*/
}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::Initialize( )
{
    SetSectionName( this->ClassName() );

    fInputHitsEvent = new TRestHitsEvent();
    fOutputHitsEvent = new TRestHitsEvent();

    fOutputEvent = fOutputHitsEvent;
    fInputEvent  = fInputHitsEvent;
}

void TRestHitsKalmanFilterProcess::LoadConfig( std::string cfgFilename, std::string name )
{

    if( LoadConfigFromFile( cfgFilename, name ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::InitProcess()
{
    //cout << __PRETTY_FUNCTION__ << endl;
}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::BeginOfEventProcess() 
{
    fOutputHitsEvent->Initialize(); 
}

//this is to get the electronVelocity
//Para E (units keV)
double TRestHitsKalmanFilterProcess::GetElectronVelocity(double E)
{
    double ElectronMassE = 0.5109989; //Mev
    double ratio = ElectronMassE *1000/ E;
    double beta = sqrt(1 - ratio * ratio);

    return beta;    
}
//this need us to check ,I think there exist some error
//here is from kalman result to get teh speed rotation
//here we need change E unit kev
//from the last state speed we get the angle
double TRestHitsKalmanFilterProcess::GetVelocityOffsetTheta02(double dZ, SVector6 XThcoli, Double_t E)
{
    double beta = GetElectronVelocity(E); //GetElectronVelocity
    double Rho = 1.662E-3;
    double ElectronMassE = 0.5109989; //Mev 
    double X0 = 19.55;
    double Pressure = 10.0;
    double X01 = 0.1 * X0 / (Pressure * Rho); // x01 is mm ? we need to check it!
    double l = sqrt(XThcoli(3)*XThcoli(3)+XThcoli(4)*XThcoli(4)+XThcoli(5)*XThcoli(5))*dZ / abs(XThcoli(5)); //mm
    double th02 = 13.6 * 13.6 / (beta *beta* (E * E * 1e-6 - ElectronMassE * ElectronMassE)) * l / X01;
    
    cout <<"dZ     "<< dZ << endl;
    cout <<"l     "<< l << endl;
    cout <<"th02     "<< th02 << endl;
    return th02;
}

//______________________________________________________________________________
TRestEvent* TRestHitsKalmanFilterProcess::ProcessEvent( TRestEvent *evInput )
{

    fInputHitsEvent = (TRestHitsEvent *) evInput;
    fOutputHitsEvent->SetEventInfo(fInputHitsEvent);
    // Copying the input hits event to the output hits event
    for( int h = 0; h < fInputHitsEvent->GetNumberOfHits(); h++ )
    {
        Double_t x = fInputHitsEvent->GetX( h );
        Double_t y = fInputHitsEvent->GetY( h );
        Double_t z = fInputHitsEvent->GetZ( h );
        Double_t en = fInputHitsEvent->GetEnergy( h );
        fOutputHitsEvent->AddHit( x, y, z, en );
    }

    //cout<< "HitsKalmanFilter_compliation_pass"<<endl;
    //1.Get the MS speed
    TRestHits *hits = fOutputHitsEvent->GetHits();
    //Get the first Hit
    Int_t Nhits=hits->GetNumberOfHits();

   
    double * MsVx=new double[Nhits];
    double * MsVy=new double[Nhits];
    double * MsVz=new double[Nhits];
    double * Energy=new double[Nhits];
    double * Speed=new double[Nhits]; // m/s
    
    //get energy and MS speed in each point
    //Double_t total_energy= fInputHitsEvent->GetEnergy(); //kev
    Double_t total_energy= 2500 + 510.9989; //kev
    //here we need the initial energy which is kinetic energy add potential energy
    //rather than deposit energy
    double InitialSpeed=lightSpeed*GetElectronVelocity(total_energy); // m/s
    Energy[0]=total_energy;
    Speed[0]=InitialSpeed;
    for( int h = 1; h <hits->GetNumberOfHits() ; h++ )
    {   //Energy[h-1] replace the true energy in h-1, 
        //hits->GetEnergy( h-1 ) is the deposit energy in h-1
        Energy[h]=Energy[h-1]-hits->GetEnergy( h-1 ); //kev
        Speed[h]=lightSpeed*GetElectronVelocity(Energy[h]);  // m/s
    } 
   
    //get 3 coordinate speed in each point
    //In the last point there no 3 coordinate speed because there no direction
    for( int h = 1; h <hits->GetNumberOfHits() ; h++ )
    {
        
       MsVx[h-1]= Speed[h-1]*( ( hits->GetX( h )-hits->GetX( h-1  ) ) / (hits->GetDistance(h,h-1)) )  ;
       MsVy[h-1]= Speed[h-1]*( ( hits->GetY( h )-hits->GetY( h-1  ) ) / (hits->GetDistance(h,h-1)) )  ;
       MsVz[h-1]= Speed[h-1]*( ( hits->GetZ( h )-hits->GetZ( h-1  ) ) / (hits->GetDistance(h,h-1)) )  ;
        //here we use the normalization speed
       MsVx[h-1]=MsVx[h-1]/Speed[h-1];
       MsVy[h-1]=MsVy[h-1]/Speed[h-1];
       MsVz[h-1]=MsVz[h-1]/Speed[h-1];
       //this is for test output
      if (fOutputHitsEvent->GetID()==0 )
        {
            cout <<"h:  "<< h \
            << "  x: "<<hits->GetX( h-1 )\
            << "  y: "<<hits->GetY( h-1 )\
            << "  z: "<<hits->GetZ( h-1 )\
            << "  MsVx["<<h-1<<"]:    "<<MsVx[h-1]\
            << "  MsVy["<<h-1<<"]:    " <<MsVy[h-1] \
            << "  MsVz["<<h-1<<"]:    "<<MsVz[h-1] \
            << "  Energy["<<h-1<<"]:    "<<Energy[h-1] \
            << "  Speed["<<h-1<<"]:    "<<Speed[h-1]<< endl;
            //cout <<"h:  "<< h<<"  x:  "<< x <<"  y:  "<< y <<"  z:  "<< z <<"  en:  "<<en<< "  hits->fX[n]:"<<  hits->fX[n] <<"n:   "<<n<< endl;
        }

    }


    //now we begin to use kalman
    SVector3 XMscoli;
    SVector6 Expectation;
    SVector6 *ExpectationRecord = new SVector6[Nhits];
    SVector6 Kalmanresult;
    SVector6 SmoothResult;

    //theoric state transition matrix F
    SMatrix66 F;
    SMatrix36 measureMatrix;
    measureMatrix(0,0) = 1.0;
    measureMatrix(1,1) = 1.0;
    measureMatrix(2,2) = 1.0;
    SMatrix66 I;
    I(0,0) = 1.0;
    I(1,1) = 1.0;
    I(2,2) = 1.0;
    I(3,3) = 1.0;
    I(4,4) = 1.0;
    I(5,5) = 1.0;

    SMatrix66 Priori;
    SMatrix66 *PrioriRecord = new SMatrix66[Nhits];
    SMatrix66 Posteriori;
    SMatrix66 *PosterioriRecord = new SMatrix66[Nhits];
    SMatrix63 KalmanGain;
    SMatrix66 SmoothGain;
    SMatrix66 NoiseThMatrix;
    SMatrix33 NoiseMsMatrix;
    SMatrix33 KIn;

    Posteriori(0,0)=1.0;
    Posteriori(1,1)=1.0;
    Posteriori(2,2)=1.0;
    Posteriori(3,3)=1.0;
    Posteriori(4,4)=1.0;
    Posteriori(5,5)=1.0;

    PosterioriRecord[0] = Posteriori;
    
    //this is for save the kalman resault.
    vector<Double_t> KFx(Nhits);
    vector<Double_t> KFy(Nhits);
    vector<Double_t> KFz(Nhits);
    vector<Double_t> KFvx(Nhits);
    vector<Double_t> KFvy(Nhits);
    vector<Double_t> KFvz(Nhits);



    // double * XMscoli=new double[Nhits];
    //double * observalueExpectationX=new double[Nhits];
    //double * observalueSpeedxX=new double[Nhits];
    //double * observalueXMscoliX=new double[Nhits];
    TRandom3 gaus;

    Double_t *xThary = new Double_t[Nhits];
    Double_t *yThary = new Double_t[Nhits];
    Double_t *zThary = new Double_t[Nhits];

    Double_t *xMsary = new Double_t[Nhits];
    Double_t *yMsary = new Double_t[Nhits];
    Double_t *zMsary = new Double_t[Nhits];

    gaus.SetSeed((unsigned) time(0));
    for (int h = 0; h < Nhits; h++)
    {
        xThary[h] = hits->GetX( h );
        yThary[h] = hits->GetY( h );
        zThary[h] = hits->GetZ( h );

        xMsary[h] = hits->GetX( h )+gaus.Gaus(0.0, MeasureError);
        yMsary[h] = hits->GetY( h )+gaus.Gaus(0.0, MeasureError);
        zMsary[h] = hits->GetZ( h )+gaus.Gaus(0.0, MeasureError);
    }

    KFx[0] = xMsary[0];
    KFy[0] = yMsary[0];
    KFz[0] = zMsary[0];
    KFvx[0] =  MsVx[0];
    KFvy[0] =  MsVy[0];
    KFvz[0] =  MsVz[0];

    Double_t thetarnd0;
    Double_t Speedx,Speedy,Speedz; // m/s
    double total_speed , dZ;
    for(Int_t n2 = 1; n2 < Nhits; n2++){
        Speedx = MsVx[n2 - 1];
        Speedy = MsVy[n2 - 1];
        Speedz = MsVz[n2 - 1];
        total_speed = sqrt(Speedx*Speedx+Speedy*Speedy+Speedz*Speedz);
       // cout<<"n2 is: "<< n2<<" " <<"Speedx is "<< Speedx << " "<<" Speedy is "<< Speedy <<" "<<"Speedz is "<< Speedz <<" "<<" total_speed is"<< total_speed<<endl ;
        Kalmanresult[0] = KFx[n2 - 1];
        Kalmanresult[1] = KFy[n2 - 1];
        Kalmanresult[2] = KFz[n2 - 1];
        Kalmanresult[3] = KFvx[n2 - 1];
        Kalmanresult[4] = KFvy[n2 - 1];
        Kalmanresult[5] = KFvz[n2 - 1];
        dZ = zThary[n2]-zThary[n2-1];
        F(0,0) = 1.0;
        F(0,3) = dZ / Speedz;
        F(1,1) = 1.0;
        F(1,4) = dZ / Speedz;
        F(2,2) = 1.0;
        F(2,5) = dZ / Speedz;
        F(3,3) = 1.0;
        F(4,4) = 1.0;
        F(5,5) = 1.0;
        
        Expectation = F * Kalmanresult;

        ExpectationRecord[n2] = Expectation;
/*
        if ( (fOutputHitsEvent->GetID()==0) )
        {   cout<<"this is lu test------:"<<endl;
            cout<<"n2 is :"<< n2 <<endl;

            cout<<"Expectation is :"<<endl;
            cout<< Expectation <<endl;

            cout<<"Kalmanresult is :"<<endl;
            cout<< Kalmanresult <<endl;
            //cout<<"Speedz is"<< Speedz <<endl;
            //cout<<"Expectation[0] is"<< Expectation[0] <<endl;

        }  
        */
        //ExpectationRecord[n2] = Expectation;

        NoiseThMatrix(0,0) = 0.0;
        NoiseThMatrix(1,1) = 0.0;
        NoiseThMatrix(2,2) = 0.0;
        cout <<"THE Kalmanresult is     " << endl;
        cout <<Kalmanresult << endl;
        thetarnd0 = sqrt(GetVelocityOffsetTheta02(dZ, Kalmanresult, Energy[n2-1]));
        cout <<"THE thetarnd0 is     "<< thetarnd0 << endl;
/*
	Double_t dx = xThary[n2] - xThary[n2 - 1];
	Double_t dy = yThary[n2] - yThary[n2 - 1];
	Double_t dz = zThary[n2] - zThary[n2 - 1];
	Double_t zhenshicos = (xThary[n2] * dx + yThary[n2] * dy + zThary[n2] * dy) / 
				(sqrt(dx *dx + dy * dy + dz * dz) * sqrt(xThary[n2] * xThary[n2] + yThary[n2] * yThary[n2] + zThary[n2] * zThary[n2]));
	Double_t zhenshi = acos(zhenshicos);
       // cout<<"THE RECURE n2 is     "<< n2 << endl;
       // cout<<"THE energe E is  "<< E << endl;
        cout <<"THE thetarnd0 is     "<< thetarnd0 << endl;
	cout <<"THE thetarnd is     "<< zhenshi << endl;
 */
       /* 
        NoiseThMatrix(3,3) = 2.0 * sin(thetarnd0/2) * sqrt(0.5 * (pow(Speedy, 2) + pow(Speedz, 2)));
        NoiseThMatrix(4,4) = 2.0 * sin(thetarnd0/2) * sqrt(0.5 * (pow(Speedx, 2) + pow(Speedz, 2)));
        NoiseThMatrix(5,5) = 2.0 * sin(thetarnd0/2) * sqrt(0.5 * (pow(Speedx, 2) + pow(Speedy, 2)));
*/
	    NoiseThMatrix(3,3) = 2.0 * sin(thetarnd0/2) * Speedx;
        NoiseThMatrix(4,4) = 2.0 * sin(thetarnd0/2) * Speedy;
        NoiseThMatrix(5,5) = 2.0 * sin(thetarnd0/2) * Speedz;

        Priori = F * Posteriori*ROOT::Math::Transpose(F) + NoiseThMatrix * NoiseThMatrix;
       /* 
        if ( (fOutputHitsEvent->GetID()==0) ){
            cout<<"thetarnd0 is :"<<endl;
            cout<< thetarnd0 <<endl;

            cout<<"sin(thetarnd0/2) is :"<<endl;
            cout<< (sin(thetarnd0/2) )<<endl;
            
            cout<<"F is :"<<endl;
            cout<< F <<endl;
            
            cout<<"NoiseThMatrix is :"<<endl;
            cout<< NoiseThMatrix <<endl;

            cout<<"Priori is :"<<endl;
            cout<< Priori <<endl;

        }
        */ 
        PrioriRecord[n2] = Priori;

        NoiseMsMatrix(0,0) = MeasureError;
        NoiseMsMatrix(1,1) = MeasureError;
        NoiseMsMatrix(2,2) = MeasureError;

        KIn = measureMatrix * Priori * ROOT::Math::Transpose(measureMatrix) + NoiseMsMatrix * NoiseMsMatrix;

        KIn.Invert();

        KalmanGain = Priori * ROOT::Math::Transpose(measureMatrix) * KIn;
        //here we need to use the theory + noise
        XMscoli[0] = xMsary[n2];
        XMscoli[1] = yMsary[n2];
        XMscoli[2] = zMsary[n2];
/*
         if ( (fOutputHitsEvent->GetID()==0) ){
            
            cout<<"KalmanGain is :"<<endl;
            cout<< KalmanGain <<endl;

        }  
  */      Kalmanresult = Expectation + KalmanGain * (XMscoli - measureMatrix * Expectation);
      /*
        if ( (n2==1 ) && (fOutputHitsEvent->GetID()==0) )
        {
            cout<<"this is lu test:"<<endl;
            cout<<"Expectation[0] is "<< Expectation[0] <<endl;
            cout<<"KalmanGain(0,0) is "<< KalmanGain(0,0) <<endl;
            cout<<"KalmanGain(0,1) is "<< KalmanGain(0,1) <<endl;
            cout<<"KalmanGain(0,2) is "<< KalmanGain(0,2) <<endl;
            cout<<"XMscoli[0] is"     << XMscoli[0] <<endl;
        }  

        if ( (n2==2 ) && (fOutputHitsEvent->GetID()==0) )
        {
            cout<<"this is lu test:"<<endl;
            cout<<"Expectation[0] is "<< Expectation[0] <<endl;
            cout<<"KalmanGain(0,0) is "<< KalmanGain(0,0) <<endl;
            cout<<"KalmanGain(0,1) is "<< KalmanGain(0,1) <<endl;
            cout<<"KalmanGain(0,2) is "<< KalmanGain(0,2) <<endl;
            cout<<"XMscoli[0] is"     << XMscoli[0] <<endl;
        } 
        */
        Posteriori = (I - KalmanGain * measureMatrix) * Priori;
        
        PosterioriRecord[n2] = Posteriori;
	
	//Add renormalization!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	double speedtotal = sqrt(Kalmanresult[3] * Kalmanresult[3] + Kalmanresult[4] * Kalmanresult[4] + Kalmanresult[5] * Kalmanresult[5]);
        KFx[n2] = Kalmanresult[0];
        KFy[n2] = Kalmanresult[1];
        KFz[n2] = Kalmanresult[2];
        KFvx[n2] = Kalmanresult[3];
        KFvy[n2] = Kalmanresult[4];
        KFvz[n2] = Kalmanresult[5];
       // KFvx[n2] = Kalmanresult[3] / speedtotal;
       // KFvy[n2] = Kalmanresult[4] / speedtotal;
       // KFvz[n2] = Kalmanresult[5] / speedtotal;
    }

    //now we begin KF Smooth
  /*  
    SmoothResult[0] = KFx[Nhits - 1];
    SmoothResult[1] = KFy[Nhits - 1];
    SmoothResult[2] = KFz[Nhits - 1];
    SmoothResult[3] = KFvx[Nhits - 1];
    SmoothResult[4] = KFvy[Nhits - 1];
    SmoothResult[5] = KFvz[Nhits - 1];
 if(fOutputHitsEvent->GetID()==0) {cout<<"the initSmooth Result Vz is "<<SmoothResult[5]<<endl;}
 
    for(Int_t n3 = Nhits - 2; n3 > -1; n3--)
    {
        Kalmanresult[0] = KFx[n3];
        Kalmanresult[1] = KFy[n3];
        Kalmanresult[2] = KFz[n3];
        Kalmanresult[3] = KFvx[n3];
        Kalmanresult[4] = KFvy[n3];
        Kalmanresult[5] = KFvz[n3];

        PrioriRecord[n3 + 1].Invert();

        //t1 -> GetEntry(n3);
        F(0,0) = 1.0;
        F(0,3) = dZ / KFvz[n3];
        F(1,1) = 1.0;
        F(1,4) = dZ / KFvz[n3];
        F(2,2) = 1.0;
        F(2,5) = dZ;
        F(3,3) = 1.0;
        F(4,4) = 1.0;
        F(5,5) = 1.0;

        SmoothGain = PosterioriRecord[n3] * ROOT::Math::Transpose(F) * PrioriRecord[n3 + 1];

         if(fOutputHitsEvent->GetID()==0)
            {    SVector6 interresult;
                interresult = SmoothResult - ExpectationRecord[n3 + 1] ;
                cout<<"SmoothGain is :"<<SmoothGain<<endl;
                cout<<"interresult is :"<<interresult<<endl;
                cout<<"PosterioriRecord["<<n3<<"] :"<<PosterioriRecord[n3]<<endl;
                cout<<"PrioriRecord["<<n3 + 1<<"] :"<<PrioriRecord[n3 + 1]<<endl;
        }

        SmoothResult = Kalmanresult + SmoothGain * (SmoothResult - ExpectationRecord[n3 + 1]);

    /*     if(fOutputHitsEvent->GetID()==0)
            {  
                
                cout<<"n3 is :"<<n3<<endl;
                cout<<"Kalmanresult is :"<<Kalmanresult<<endl;
                cout<<"SmoothResult is :"<<SmoothResult<<endl;
                cout<<"ExpectationR is :"<<ExpectationRecord[n3 + 1]<<endl;
                
                cout<<"SmoothGain is :"<<SmoothGain<<endl;
                
                
                cout<<"SmoothResult[2] is :"<<SmoothResult[2]<<endl;
                //
        }




	//?????????????????????????????
       //last zhushi 
        KFx[n3] = SmoothResult[0];
        KFy[n3] = SmoothResult[1];
        KFz[n3] = SmoothResult[2];
        KFvx[n3] = SmoothResult[3];
        KFvy[n3] = SmoothResult[4];
        KFvz[n3] = SmoothResult[5];
    }
    */
    //there still have kalman speed information But have not been saved.
    for (int h = 0; h < Nhits; h++)
    {
        hits->fX[h]=KFx[h];
        hits->fY[h]=KFy[h];
        hits->fZ[h]=KFz[h];
            
    }
  /*  for (int h = 0; h < Nhits; h++){

      if (isnan( KFz[h] ) == 1 )
        {
            return NULL;
        }  
            
    }
    */
    //this is for test output 
    if (fOutputHitsEvent->GetID()==0 )
    {
       cout << "output event" << endl;
       cout << "+++++++++++++++++" << endl;
       fOutputHitsEvent->PrintEvent();
       cout << "+++++++++++++++++" << endl;
    }

    //draw 
    //if (fOutputHitsEvent->GetID()==0 )
     {
    Double_t *xKfary = new Double_t[Nhits];
    Double_t *yKfary = new Double_t[Nhits];
    Double_t *zKfary = new Double_t[Nhits];

    ofstream outfileThMs;
    string str1="___ThMsposition.txt";
    outfileThMs.open(to_string(fOutputHitsEvent->GetID())+str1);
    outfileThMs << Nhits <<endl;
    for (int h = 0; h < Nhits; h++)
    {
      
        xKfary[h]= KFx[h];
        yKfary[h]= KFy[h];
        zKfary[h]= KFz[h];

    outfileThMs << Energy[h] << " " \
                << xThary[h] << " " << yThary[h] << " " << zThary[h] << " " \
                << xMsary[h] << " " << yMsary[h] << " " << zMsary[h] << " " \
                << xKfary[h] << " " << yKfary[h] << " " << zKfary[h] <<endl;
    
    }
    outfileThMs.close();
 /*   
    TCanvas *c1=new TCanvas("cXY","cXY",800,800);
 
    TGraph *gr11 = new TGraph (Nhits-3, xMsary, yMsary);
    gr11 -> GetXaxis()->SetTitle("X [mm]");
    gr11 -> GetYaxis()->SetTitle("Y [mm]");
    gr11 -> GetYaxis()->SetTitleOffset(1.4);
    gr11 -> SetTitle("XY Projection");
    gr11 -> SetLineWidth(3);
    gr11 -> SetLineColor(4);
    gr11 -> Draw("AC");
    /*
    TGraph *gr12 = new TGraph (Nhits-10, xMsary, yMsary);
    gr12 -> GetXaxis()->SetTitle("X [mm]");
    gr12 -> GetYaxis()->SetTitle("Y [mm]");
    //gr12 -> SetMarkerStyle(7);
    //gr12 -> SetMarkerSize(3);
    //gr12 -> SetMarkerColor(1);
    gr12 -> SetLineWidth(3);
    gr12 -> SetLineColor(16);
    //gr12 -> Draw("P""same");
    gr12 -> Draw("PC");
    */
 /*   TGraph *gr13 = new TGraph (Nhits-3, xKfary, yKfary);
    gr13 -> SetMarkerStyle(7);
    gr13 -> SetMarkerSize(3);
    gr13 -> SetMarkerColor(2);
    gr13 -> SetLineWidth(3);
    gr13 -> SetLineColor(7);
    //gr13 -> Draw("P""same");
    gr13 -> Draw("PC");
/*
    TLegend *legend1 = new TLegend(0.1,0.75,0.35,0.9); 
    //legend1 -> AddEntry(gr11,"Simulation Data","l");
    legend1 -> AddEntry(gr11,"Measure Data","p");
    legend1 -> AddEntry(gr13,"KF Result","p");
    legend1 -> Draw("same");
   //cout<<"show the pic"<<endl;
    //getchar();
    c1->Print("Event0_kf.png");
*/
    delete[] xThary;
    delete[] yThary;
    delete[] zThary;

    delete[] xMsary;
    delete[] yMsary;
    delete[] zMsary;

    delete[] xKfary;
    delete[] yKfary;
    delete[] zKfary;


    //getchar();
    }


    return fOutputHitsEvent;
}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::EndOfEventProcess() 
{

}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::EndProcess()
{

}

//______________________________________________________________________________
void TRestHitsKalmanFilterProcess::InitFromConfigFile( )
{
 /*   fStartingDistance = GetDblParameterWithUnits(  "startingDistance" );
    fMinimumDistance  = GetDblParameterWithUnits( "minimumDistance" );
    fDistanceFactor   = StringToDouble( GetParameter( "distanceStepFactor" ) );
    fMaxNodes         = StringToDouble( GetParameter( "maxNodes" ) );
*/
}

