///______________________________________________________________________________
///______________________________________________________________________________
///             
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestG4toHitsGetPhotoElectricProcess.cxx
///
///
///             Simple process to convert a TRestG4Event class into a 
///    		    TRestHitsEvent, that is, we just "extract" the hits information
///             Date : oct/2016
///             Author : I. G. Irastorza
///
///_______________________________________________________________________________


#include "TRestG4toHitsGetPhotoElectricProcess.h"
using namespace std;


ClassImp(TRestG4toHitsGetPhotoElectricProcess)
    //______________________________________________________________________________
TRestG4toHitsGetPhotoElectricProcess::TRestG4toHitsGetPhotoElectricProcess()
{
    Initialize();
}

//______________________________________________________________________________
TRestG4toHitsGetPhotoElectricProcess::TRestG4toHitsGetPhotoElectricProcess( char *cfgFileName )
{
    Initialize();

    if( LoadConfigFromFile( cfgFileName ) ) LoadDefaultConfig( );
}

//______________________________________________________________________________
TRestG4toHitsGetPhotoElectricProcess::~TRestG4toHitsGetPhotoElectricProcess()
{
    delete fG4Event;
    delete fHitsEvent;
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::LoadDefaultConfig()
{
    SetTitle( "Default config" );

    cout << "G4 to hits metadata not found. Loading default values" << endl;
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::Initialize()
{
    SetSectionName( this->ClassName() );

    fG4Event = new TRestG4Event();
    fHitsEvent = new TRestHitsEvent();

    fOutputEvent = fHitsEvent;
    fInputEvent = fG4Event;
}

void TRestG4toHitsGetPhotoElectricProcess::LoadConfig( std::string cfgFilename, std::string name )
{
    if( LoadConfigFromFile( cfgFilename, name ) ) LoadDefaultConfig( );
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::InitProcess()
{
    //    TRestEventProcess::ReadObservables();

    fG4Metadata = (TRestG4Metadata *) GetGeant4Metadata( );

    for( unsigned int n = 0; n < fVolumeSelection.size(); n++ )
    {
        if( fG4Metadata->GetActiveVolumeID( fVolumeSelection[n] ) >= 0 ) 
            fVolumeId.push_back( fG4Metadata->GetActiveVolumeID( fVolumeSelection[n] ) ); 
        else if( GetVerboseLevel() >= REST_Warning )
            cout << "TRestG4toHitsGetPhotoElectricProcess. volume name : " << fVolumeSelection[n] << " not found and will not be added." << endl;
    }
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::BeginOfEventProcess() 
{
    fHitsEvent->Initialize();
}

//______________________________________________________________________________
TRestEvent* TRestG4toHitsGetPhotoElectricProcess::ProcessEvent( TRestEvent *evInput )
{
    fG4Event = (TRestG4Event *) evInput;
    if (fG4Event->GetTotalDepositedEnergy()<2300)
    {
        return NULL;
    }
    if ( ! ( fG4Event->isPhotoElectric( ) ) )
    {
        return NULL;
    }
    
    fHitsEvent->SetRunOrigin( fG4Event->GetRunOrigin() );
    fHitsEvent->SetSubRunOrigin( fG4Event->GetSubRunOrigin() );
    fHitsEvent->SetID( fG4Event->GetID() );
    fHitsEvent->SetSubID( fG4Event->GetSubID() );
    fHitsEvent->SetSubEventTag( fG4Event->GetSubEventTag() );
    fHitsEvent->SetTimeStamp( fG4Event->GetTimeStamp() );
    fHitsEvent->SetState( fG4Event->isOk() );

    Int_t i,j;
    Double_t x,y,z,E;
    Int_t TrackIsPhotoElectricNum = 0;
    Int_t TrackIsComptonNum = 0;
    Int_t TrackIsRadiactiveDecay = 0;

    for ( i = 0; i < fG4Event->GetNumberOfTracks(); i++ )
    { 
        /*
        if( fG4Event->GetTrack( i )->isPhotoElectric( ) ) TrackIsPhotoElectricNum+=1; 
        if( fG4Event->GetTrack( i )->isCompton( ) ) TrackIsComptonNum+=1; 
        if( fG4Event->GetTrack( i )->isRadiactiveDecay( ) ) TrackIsRadiactiveDecay+=1; 
        if( !( fG4Event->GetTrack( i )->isPhotoElectric( ) )  ) continue;
        */

        //if( fG4Event->GetTrack( i )->isCompton( ) ) continue;        
        for ( j = 0; j < fG4Event->GetTrack(i)->GetNumberOfHits(); j++ )
        {
           // read x,y,z and E of every hit in the G4 event
            x = fG4Event->GetTrack(i)->GetHits()->fX[j];
            y = fG4Event->GetTrack(i)->GetHits()->fY[j];
            z = fG4Event->GetTrack(i)->GetHits()->fZ[j];
            E = fG4Event->GetTrack(i)->GetHits()->fEnergy[j];

            Bool_t addHit = true;
            if( fVolumeId.size() > 0 )
            {
                addHit = false;
                for( unsigned int n = 0; n < fVolumeId.size(); n++ )
                    if( fG4Event->GetTrack(i)->GetHits()->GetVolumeId(j) == fVolumeId[n] )
                        addHit = true;
            }

            // and write them in the output hits event:
            if( addHit && E > 0 ) fHitsEvent->AddHit (x, y, z, E);
        }
    }
    


    /*
    if (fG4Event->GetID()==0)
    {
      cout << "Event ID:   "               <<fG4Event->GetID()<<endl;
      cout << "trackNUM:   "<<fG4Event->GetNumberOfTracks()<<endl;
      cout << "TrackIsPhotoElectricNum:   "<<TrackIsPhotoElectricNum<<endl;
      cout << "TrackIsComptonNum:   "<<TrackIsComptonNum<<endl;
      cout << "TrackIsRadiactiveDecay:   "<<TrackIsRadiactiveDecay<<endl;
    }
    */
    if( this->GetVerboseLevel() >= REST_Debug ) 
    {
        cout << "TRestG4toHitsGetPhotoElectricProcess. Hits added : " << fHitsEvent->GetNumberOfHits() << endl;
        cout << "TRestG4toHitsGetPhotoElectricProcess. Hits total energy : " << fHitsEvent->GetEnergy() << endl;
    }

    return fHitsEvent;
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::EndOfEventProcess() 
{

}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::EndProcess()
{
    // Function to be executed once at the end of the process 
    // (after all events have been processed)

    //Start by calling the EndProcess function of the abstract class. 
    //Comment this if you don't want it.
    //TRestEventProcess::EndProcess();
}

//______________________________________________________________________________
void TRestG4toHitsGetPhotoElectricProcess::InitFromConfigFile( )
{
    size_t position = 0;
    string addVolumeDefinition;

    while( ( addVolumeDefinition = GetKEYDefinition( "addVolume", position ) ) != "" )
        fVolumeSelection.push_back( GetFieldValue( "name", addVolumeDefinition ) );
}
