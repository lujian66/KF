///______________________________________________________________________________
///______________________________________________________________________________
///______________________________________________________________________________
///             
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestHitsMergeProcess.cxx
///
///             Jan 2016:   First concept (Javier Galan)
//
///_______________________________________________________________________________

#include "TRestHitsMergeProcess.h"
using namespace std;

ClassImp(TRestHitsMergeProcess)
//______________________________________________________________________________
TRestHitsMergeProcess::TRestHitsMergeProcess( )
{
    Initialize();
}

//______________________________________________________________________________
TRestHitsMergeProcess::TRestHitsMergeProcess( char *cfgFileName )
{
    Initialize();

    if( LoadConfigFromFile( cfgFileName ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
TRestHitsMergeProcess::~TRestHitsMergeProcess( )
{
    delete fInputHitsEvent;
    delete fOutputHitsEvent;
}

void TRestHitsMergeProcess::LoadDefaultConfig( )
{
    SetName( "hitsReductionProcess" );
    SetTitle( "Default config" );

    fStartingDistance = 0.5;
    fMinimumDistance  = 3;
    fDistanceFactor   = 1.5;
    fMaxNodes         = 30;
}

//______________________________________________________________________________
void TRestHitsMergeProcess::Initialize( )
{
    SetSectionName( this->ClassName() );

    fInputHitsEvent = new TRestHitsEvent();
    fOutputHitsEvent = new TRestHitsEvent();

    fOutputEvent = fOutputHitsEvent;
    fInputEvent  = fInputHitsEvent;
}

void TRestHitsMergeProcess::LoadConfig( std::string cfgFilename, std::string name )
{

    if( LoadConfigFromFile( cfgFilename, name ) == -1 ) LoadDefaultConfig( );
}

//______________________________________________________________________________
void TRestHitsMergeProcess::InitProcess()
{
    cout << __PRETTY_FUNCTION__ << endl;
}

//______________________________________________________________________________
void TRestHitsMergeProcess::BeginOfEventProcess() 
{
    fOutputHitsEvent->Initialize(); 
}

//______________________________________________________________________________
TRestEvent* TRestHitsMergeProcess::ProcessEvent( TRestEvent *evInput )
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

    // Reducing the hits
    TRestHits *hits = fOutputHitsEvent->GetHits();



    Double_t distance = fStartingDistance;
    //I don't need to merge hits less than fMaxNodes, just merge the hits in 3 mm
    while( distance < fMinimumDistance  )
    {
        Bool_t merged = true;
        while( merged )
        {
            merged = false;
            for( int i = 0; i < hits->GetNumberOfHits(); i++ )
            {
                for( int j = i+1; j < hits->GetNumberOfHits(); j++ )
                {
                    if( hits->GetDistance2( i, j ) < distance * distance )
                    {
                        hits->MergeHits( i, j );
                        merged = true;
                    }
                }
            }
        }
        distance *= fDistanceFactor;
        cout<<"distance is "<<distance <<endl;
    }

    Int_t initialHits = fInputHitsEvent->GetNumberOfHits();
    Int_t finalHits = fOutputHitsEvent->GetNumberOfHits();

    cout << "TRestHitsMergeProcess : Initial number of hits : " << initialHits << endl;
    cout << "TRestHitsMergeProcess : Final number of hits : " << finalHits << endl;
    cout << "fStartingDistance is : " << fStartingDistance << endl;
     cout << "fMinimumDistance is : " << fMinimumDistance << endl;
      cout << "fDistanceFactor is : " << fDistanceFactor << endl;
        cout << "fMaxNodes is : " << fMaxNodes << endl;
    
    if( this->GetVerboseLevel() == REST_Debug )
    {
        cout << "TRestHitsMergeProcess : Initial number of hits : " << initialHits << endl;
        cout << "TRestHitsMergeProcess : Final number of hits : " << finalHits << endl;
    }

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
void TRestHitsMergeProcess::EndOfEventProcess() 
{

}

//______________________________________________________________________________
void TRestHitsMergeProcess::EndProcess()
{
}

//______________________________________________________________________________
void TRestHitsMergeProcess::InitFromConfigFile( )
{
    fStartingDistance = GetDblParameterWithUnits(  "startingDistance" );
    fMinimumDistance  = GetDblParameterWithUnits( "minimumDistance" );
    fDistanceFactor   = StringToDouble( GetParameter( "distanceStepFactor" ) );
    fMaxNodes         = StringToDouble( GetParameter( "maxNodes" ) );
}

