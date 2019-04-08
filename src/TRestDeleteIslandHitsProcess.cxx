///______________________________________________________________________________
///______________________________________________________________________________
///______________________________________________________________________________
///             
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestDeleteIslandHitsProcess.cxx
///
///             MAY 2019 : lujian
//
///_______________________________________________________________________________
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <float.h>
#include "TRestDeleteIslandHitsProcess.h"
//#include <fstream>
using namespace std;

ClassImp(TRestDeleteIslandHitsProcess)
//______________________________________________________________________________
TRestDeleteIslandHitsProcess::TRestDeleteIslandHitsProcess( )
{
    Initialize();
}

//______________________________________________________________________________
TRestDeleteIslandHitsProcess::TRestDeleteIslandHitsProcess( char *cfgFileName )
{
    Initialize();

    if( LoadConfigFromFile( cfgFileName ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
TRestDeleteIslandHitsProcess::~TRestDeleteIslandHitsProcess( )
{
    delete fInputHitsEvent;
    delete fOutputHitsEvent;
}

void TRestDeleteIslandHitsProcess::LoadDefaultConfig( )
{
    SetName( "DeleteIslandHits" );
    SetTitle( "Default config" );

  /*  fStartingDistance = 0.5;
    fMinimumDistance  = 3;
    fDistanceFactor   = 1.5;
    fMaxNodes         = 30;
*/
}

//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::Initialize( )
{
    SetSectionName( this->ClassName() );

    fInputHitsEvent = new TRestHitsEvent();
    fOutputHitsEvent = new TRestHitsEvent();

    fOutputEvent = fOutputHitsEvent;
    fInputEvent  = fInputHitsEvent;
}

void TRestDeleteIslandHitsProcess::LoadConfig( std::string cfgFilename, std::string name )
{

    if( LoadConfigFromFile( cfgFilename, name ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::InitProcess()
{
    //cout << __PRETTY_FUNCTION__ << endl;
}

//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::BeginOfEventProcess() 
{
    fOutputHitsEvent->Initialize(); 
}

/*
void TRestDeleteIslandHitsProcess::GetMaxZ(TRestEvent *evInput) 
{
    fOutputHitsEvent->Initialize(); 
}

void TRestDeleteIslandHitsProcess::GetMinZ(TRestEvent *evInput) 
{
    fOutputHitsEvent->Initialize(); 
}
*/
Double_t TRestDeleteIslandHitsProcess::GetMax(Double_t *ary, Int_t N)
{
    Double_t Max = -9999999.9;
    for(Int_t i = 0; i < N; i++)
    {
        if (Max < ary[i])
        {
            Max = ary[i];
        }
    }
    return Max;
}

Double_t TRestDeleteIslandHitsProcess::GetMin(Double_t *ary, Int_t N)
{
    Double_t Min = 99999999.9;
    for(Int_t i = 0; i < N; i++)
    {
        if (Min > ary[i])
        {
            Min = ary[i];
        }
    }
    return Min;
}

Double_t TRestDeleteIslandHitsProcess::GetZMinValue(Double_t z,int index, double * ZArray, Int_t ZArraySize)
{
    Double_t Min = 99999999.9;
    Double_t zGap=0;
    for(Int_t i = 0; i < ZArraySize; i++)
    { 
        if (i==index) continue;
        zGap=TMath::Abs(z-ZArray[i]);
        if (Min > zGap)
        {
            Min = zGap;
        }
    }
    return Min;
}
//there need one array and array size
/*
double TRestDeleteIslandHitsProcess::Draw1DHist(double * Array, int ArraySize ) 
{
    TCanvas *c = new TCanvas("k","khist",800,800);

    Double_t ArrayMax = GetMax(Array, ArraySize);
    Double_t ArrayMin = GetMin(Array, ArraySize);

      TH1F *h4 = new TH1F("h4", "", ArraySize, ArrayMin, ArrayMax);
      h4 -> SetTitle("khist");
      //h4 -> SetStats(0);
      h4 -> GetXaxis() -> SetTitle("kvalue");
      h4 -> GetXaxis() -> SetTitleOffset(1.5);
      h4 -> GetYaxis() -> SetTitle("num");
      h4 -> GetYaxis() -> SetTitleOffset(1.4);
      for(Int_t i = 0; i < ArraySize; i++)
      {
        h4->Fill(Array[i]);
      }
      h4 -> SetMarkerStyle(21);
      h4 -> SetMarkerSize(1);
      h4 -> SetLineWidth(2);
      h4 -> SetLineColor(1);
      h4 -> Draw();
     // double khistmean = h4->GetMean();
     // cout<<"khistmean is :"<<khistmean<<endl;

     // TSpectrum *s1 = new TSpectrum(2);
//Int_t nkfound = s1->Search(h4,2,"",0.001);
      Double_t *kpeaks;
    //  kpeaks=s1->GetPositionX();
      double peakValue=kpeaks[0];
      //getchar();
      //cout<<"this first peak is "<<kpeaks[0]<<endl;
      //cout<<"this second peak is "<<kpeaks[1]<<endl;
      cout<<"process is come to here"<<endl;
      c->Print("1.png");
      delete h4;
      delete c;
      //delete s1;
     // delete kpeaks;
      return peakValue;
}
*/
//______________________________________________________________________________
TRestEvent* TRestDeleteIslandHitsProcess::ProcessEvent( TRestEvent *evInput )
{

    fInputHitsEvent = (TRestHitsEvent *) evInput;
    if(fInputHitsEvent->GetEnergy() < 2450) 
        return NULL;
    fOutputHitsEvent->SetEventInfo(fInputHitsEvent);
    TRestHits *hits = fInputHitsEvent->GetHits();
/*
    int totalHits = fInputHitsEvent->GetNumberOfHits();

    int ParaNum=totalHits*(totalHits-1)/2;
      double *k=new double [ParaNum];
      double *b=new double [ParaNum];
      double k0;
      double b0;

    double *ArrayX =new double [totalHits] ;
    double *ArrayY =new double [totalHits] ;
    double *ArrayZ =new double [totalHits] ;
    double *ArrayEn =new double [totalHits] ;
   for( int h = 0; h <fInputHitsEvent->GetNumberOfHits() ; h++ )
    {
        

        ArrayX [h] = fInputHitsEvent->GetX( h );
        ArrayY [h] = fInputHitsEvent->GetY( h );
        ArrayZ [h] = fInputHitsEvent->GetZ( h );
        ArrayEn[h] = fInputHitsEvent->GetEnergy( h );
        //h is index
        //ZMin=GetZMinValue(z , h , ZArray , totalHits);
       
       //if(  (z >ZMax3sigma )  || (z <ZMin3sigma ) ) continue;
       //fOutputHitsEvent->AddHit (x, y, z, en);
     }
  int fenmuequ0=0;
  int loopnum=0;
  for (int i = 0; i < totalHits; i++)
  {

    for (int j = i+1; j < totalHits; j++)
    {
      if (ArrayX[j]-ArrayX[i]==0)
      {
        fenmuequ0+=1;
        continue;
      }
      loopnum+=1;
       k0=(ArrayZ[j]-ArrayZ[i])/(ArrayX[j]-ArrayX[i]);
       b0=ArrayZ[i]-k0*ArrayX[i];
       k[loopnum-1]=k0;
       b[loopnum-1]=b0;
    }
   
  }
  Double_t kvMax = GetMax(k, ParaNum);
  Double_t kvMin = GetMin(k, ParaNum);

  TH1F *h4 = new TH1F("h4", "", ParaNum, kvMin, kvMax);
  for(Int_t i = 0; i < ParaNum; i++)
  {
    h4->Fill(k[i]);
  }
  TSpectrum *s1 = new TSpectrum(2);
  Int_t nkfound = s1->Search(h4,2,"",0.001);
  Double_t *kpeaks;
  kpeaks=s1->GetPositionX();
  cout<<"this first peak is "<<kpeaks[0]<<endl;
  cout<<"this second peak is "<<kpeaks[1]<<endl;


  Double_t bvMax = GetMax(b, ParaNum);
  Double_t bvMin = GetMin(b, ParaNum);
  TH1F *h5 = new TH1F("h5", "", 300, bvMin, bvMax);
  for(Int_t i = 0; i < ParaNum; i++)
  {
    h5->Fill(b[i]);
  }
  TSpectrum *s2 = new TSpectrum(2);
  Int_t nbfound = s2->Search(h5,2,"",0.001);
  Double_t *bpeaks;
  bpeaks=s2->GetPositionX();
  cout<<"this first peak is "<<bpeaks[0]<<endl;
  cout<<"this second peak is "<<bpeaks[1]<<endl;
    double *ArrayX =new double [totalHits] ;
    double *ArrayY =new double [totalHits] ;
    double *ArrayZ =new double [totalHits] ;
    double *ArrayEn =new double [totalHits] ;
   for( int h = 0; h <fInputHitsEvent->GetNumberOfHits() ; h++ )
    {
        

        ArrayX [h] = fInputHitsEvent->GetX( h );
        ArrayY [h] = fInputHitsEvent->GetY( h );
        ArrayZ [h] = fInputHitsEvent->GetZ( h );
        ArrayEn[h] = fInputHitsEvent->GetEnergy( h );
        //h is index
        //ZMin=GetZMinValue(z , h , ZArray , totalHits);
       
       //if(  (z >ZMax3sigma )  || (z <ZMin3sigma ) ) continue;
       //fOutputHitsEvent->AddHit (x, y, z, en);
     }
//sort
     double tempX,tempY,tempZ ,tempEn;
     for (int i = 0; i < totalHits-1; i++)
     {
        for (int j = 0; j <totalHits-1-i ; j++)
        {
          if (ArrayZ [j] > ArrayZ [j+1] )
          {
             tempZ = ArrayZ [j+1];
             tempY = ArrayY [j+1];
             tempX = ArrayX [j+1];
             tempEn = ArrayEn[j+1];

             ArrayZ [j+1] = ArrayZ [j];
             ArrayY [j+1] = ArrayY [j];
             ArrayX [j+1] = ArrayX [j];
             ArrayEn[j+1] = ArrayEn[j];
             
             ArrayZ [j] = tempZ ;
             ArrayY [j] = tempY ;
             ArrayX [j] = tempX ;
             ArrayEn[j] = tempEn;
          }
        }
     }

    //calculate the Gap beteween Z
    double *zGAP = new double[totalHits-1] ;
    for (int i = 0; i < totalHits-1; i++)
    {
      zGAP[i] = ArrayZ [i+1] - ArrayZ [i];
     if (fInputHitsEvent->GetID()==0) cout<<"zGAP["<<i<<"]"<<zGAP[i]<<endl;
    }
    
     for (int j = 0; j < totalHits-1; j++)
     {
       zGAP=ArrayZ [j+1] - ArrayZ [j];
         if (zGAP > 6 )
        {
            return NULL;
        }
     }

    for( int h = 0; h <fInputHitsEvent->GetNumberOfHits() ; h++ )
    {

       fOutputHitsEvent->AddHit (ArrayX [h], ArrayY [h], ArrayZ [h], ArrayEn[h]);
     }
     */
//there is to use the 3sigma to delete the outlier, 
//but there still have other way to do this things
    //TRestHits *hits = fInputHitsEvent->GetHits();

    int totalHits = fInputHitsEvent->GetNumberOfHits();
    double * ZArray=new double [totalHits];
    double ZMean, Zsigma, ZtotalGap2=0, Ztotal=0;
    for (int i = 0; i < totalHits; i++){ 
        ZArray[i]=hits->fZ[i];
        Ztotal+=ZArray[i];
    }


    ZMean=Ztotal/totalHits;

     for (int i = 0; i < totalHits; i++){ 
        ZtotalGap2+=(ZArray[i] - ZMean) * (ZArray[i] - ZMean);
    }

    Zsigma=TMath::Sqrt(ZtotalGap2/totalHits);

    if (fInputHitsEvent->GetID()==0)
    {
       cout << "output event" << endl;
       cout << "fInputHitsEvent->GetID() is" <<fInputHitsEvent->GetID() << endl;
       cout << "+++++++++++++++++" << endl;
       for (int i = 0; i < 5; i++)
       {
           cout<<"ZArray["<<i<<"] is"<<ZArray[i]<<endl;
       }
       cout << "+++++++++++++++++" << endl;
    }
    double threshold=6 ; //mm it means the Z Gap more than 6 than we though it belongs to island 
    double ZMin3sigma=ZMean-Zsigma*3;
    double ZMax3sigma=ZMean+Zsigma*3;
    for( int h = 0; h <fInputHitsEvent->GetNumberOfHits() ; h++ )
    {
        

        Double_t x = fInputHitsEvent->GetX( h );
        Double_t y = fInputHitsEvent->GetY( h );
        Double_t z = fInputHitsEvent->GetZ( h );
        Double_t en = fInputHitsEvent->GetEnergy( h );
        //h is index
        //ZMin=GetZMinValue(z , h , ZArray , totalHits);
       
       if(  (z >ZMax3sigma )  || (z <ZMin3sigma ) ) continue;
       fOutputHitsEvent->AddHit (x, y, z, en);
     }

    if (fInputHitsEvent->GetID()==0)
    {
       cout << "output event" << endl;
       cout << "fInputHitsEvent->GetID() is" <<fInputHitsEvent->GetID() << endl;
       cout << "+++++++++++++++++" << endl;
       //for (int i = 0; i < totalHits; i++)
       //{
       //    cout<<"h:""ZArray["<<i<<"] is"<<ZArray[i]<<endl;
       //}
       fOutputHitsEvent->PrintEvent();
       cout << "+++++++++++++++++" << endl;
    }
    Int_t initialHits = fInputHitsEvent->GetNumberOfHits();
    Int_t finalHits = fOutputHitsEvent->GetNumberOfHits();

    cout << "TRestHitsdeleteProcess : Initial number of hits : " << initialHits << endl;
    cout << "TRestHitsdeleteProcess : Final number of hits : " << finalHits << endl;
  //  delete [] ArrayX;
  //  delete [] ArrayY;
  //  delete [] ArrayZ;
  //  delete [] ArrayEn;



return fOutputHitsEvent;
}




//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::EndOfEventProcess() 
{

}

//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::EndProcess()
{
}

//______________________________________________________________________________
void TRestDeleteIslandHitsProcess::InitFromConfigFile( )
{
 /*   fStartingDistance = GetDblParameterWithUnits(  "startingDistance" );
    fMinimumDistance  = GetDblParameterWithUnits( "minimumDistance" );
    fDistanceFactor   = StringToDouble( GetParameter( "distanceStepFactor" ) );
    fMaxNodes         = StringToDouble( GetParameter( "maxNodes" ) );
*/
}

