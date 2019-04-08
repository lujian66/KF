///______________________________________________________________________________
///______________________________________________________________________________
///______________________________________________________________________________
///             
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestHitsPretreatmentProcess.cxx
///
///             MAY 2019 : lujian
//
///_______________________________________________________________________________
#include <float.h>
#include "TRestHitsPretreatmentProcess.h"
using namespace std;

ClassImp(TRestHitsPretreatmentProcess)
//______________________________________________________________________________
TRestHitsPretreatmentProcess::TRestHitsPretreatmentProcess( )
{
    Initialize();
}

//______________________________________________________________________________
TRestHitsPretreatmentProcess::TRestHitsPretreatmentProcess( char *cfgFileName )
{
    Initialize();

    if( LoadConfigFromFile( cfgFileName ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
TRestHitsPretreatmentProcess::~TRestHitsPretreatmentProcess( )
{
    delete fInputHitsEvent;
    delete fOutputHitsEvent;
}

void TRestHitsPretreatmentProcess::LoadDefaultConfig( )
{
    SetName( "hitsReductionProcess" );
    SetTitle( "Default config" );

  /*  fStartingDistance = 0.5;
    fMinimumDistance  = 3;
    fDistanceFactor   = 1.5;
    fMaxNodes         = 30;
*/
}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::Initialize( )
{
    SetSectionName( this->ClassName() );

    fInputHitsEvent = new TRestHitsEvent();
    fOutputHitsEvent = new TRestHitsEvent();

    fOutputEvent = fOutputHitsEvent;
    fInputEvent  = fInputHitsEvent;
}

void TRestHitsPretreatmentProcess::LoadConfig( std::string cfgFilename, std::string name )
{

    if( LoadConfigFromFile( cfgFilename, name ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::InitProcess()
{
    //cout << __PRETTY_FUNCTION__ << endl;
}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::BeginOfEventProcess() 
{
    fOutputHitsEvent->Initialize(); 
}

/*
void TRestHitsPretreatmentProcess::GetMaxZ(TRestEvent *evInput) 
{
    fOutputHitsEvent->Initialize(); 
}

void TRestHitsPretreatmentProcess::GetMinZ(TRestEvent *evInput) 
{
    fOutputHitsEvent->Initialize(); 
}
*/

//______________________________________________________________________________
TRestEvent* TRestHitsPretreatmentProcess::ProcessEvent( TRestEvent *evInput )
{

    fInputHitsEvent = (TRestHitsEvent *) evInput;
    if(fInputHitsEvent->GetEnergy() < 2000) 
        return NULL;
    fOutputHitsEvent->SetEventInfo(fInputHitsEvent);

    Double_t ZMax= -1500.0;
    Double_t ZMin=  1500.0;
    // get the Zmax and Zmin

    cout<<"ZMax:    "<<ZMax<<endl;
    cout<<"ZMin:    "<<ZMin<<endl;

    for( int h = 0; h < fInputHitsEvent->GetNumberOfHits(); h++ )
    {
        Double_t z = fInputHitsEvent->GetZ( h );
        if(ZMin>z) ZMin=z;
        if(ZMax<z) ZMax=z;
    }
    cout<<"ZMax_the_second:    "<<ZMax<<endl;
    cout<<"ZMin_the_second:    "<<ZMin<<endl;

    TRestHits *hits = fOutputHitsEvent->GetHits();
      
    Int_t Nhits;

    Nhits= (ZMax-(ZMin-1.5) )/3+1;

    cout<<"Nhits:    "<<Nhits<<endl;

    for (int i = 0; i <Nhits; i++)
    {
        hits->AddHit(0,0,0,0);
    }

    Int_t n=0; 
    for( int h = 0; h <fInputHitsEvent->GetNumberOfHits() ; h++ )
    {
        

        Double_t x = fInputHitsEvent->GetX( h );
        Double_t y = fInputHitsEvent->GetY( h );
        Double_t z = fInputHitsEvent->GetZ( h );
        Double_t en = fInputHitsEvent->GetEnergy( h );
        
        n=(z-(ZMin-1.5))/3;

        hits->fX[n]=(x*en+(hits->fX[n])*(hits->fEnergy[n]))/(en+hits->fEnergy[n]);
        hits->fY[n]=(y*en+(hits->fY[n])*(hits->fEnergy[n]))/(en+hits->fEnergy[n]);
        hits->fZ[n]=ZMin+n*3;

        hits->fEnergy[n]=hits->fEnergy[n]+en;

       //if (fOutputHitsEvent->GetID()==1 )
       // {
            //cout <<"h:  "<< h << "  hits->fEnergy[n]:"<<  hits->fEnergy[n] <<"n:   "<<n<< endl;
            cout <<"h:  "<< h<<"  x:  "<< x <<"  y:  "<< y <<"  z:  "<< z <<"  en:  "<<en<< "  hits->fX[n]:"<<  hits->fX[n] <<"n:   "<<n<< endl;
       // }
       //cout << "fZ[n]:"<<  hits->fZ[n] << endl;

    }

    Double_t total_energy=0;

    for (int i = 0; i < Nhits; i++)
    {
        total_energy+=hits->fEnergy[i];
    }

    hits->fTotEnergy=total_energy;
   
    cout << "TRestHitsPretreatmentProcess : total_energy : " << total_energy << endl;


    Int_t initialHits = fInputHitsEvent->GetNumberOfHits();
    Int_t finalHits = fOutputHitsEvent->GetNumberOfHits();

   // if( this->GetVerboseLevel() == REST_Debug )
   // {
        cout << "TRestHitsPretreatmentProcess : Initial number of hits : " << initialHits << endl;
        cout << "TRestHitsPretreatmentProcess : Final number of hits : " << finalHits << endl;
  //  }

    /*
       cout << "output event" << endl;
       cout << "+++++++++++++++++" << endl;
       fOutputHitsEvent->PrintEvent();
       cout << "+++++++++++++++++" << endl;
       getchar();
    */


    //cout << "Number output of tracks : " << fOutputTrackEvent->GetNumberOfTracks() << endl;

    return fOutputHitsEvent;
}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::EndOfEventProcess() 
{

}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::EndProcess()
{
}

//______________________________________________________________________________
void TRestHitsPretreatmentProcess::InitFromConfigFile( )
{
 /*   fStartingDistance = GetDblParameterWithUnits(  "startingDistance" );
    fMinimumDistance  = GetDblParameterWithUnits( "minimumDistance" );
    fDistanceFactor   = StringToDouble( GetParameter( "distanceStepFactor" ) );
    fMaxNodes         = StringToDouble( GetParameter( "maxNodes" ) );
*/
}

