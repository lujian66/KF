///______________________________________________________________________________
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
///             TRestHitsEvent, that is, we just "extract" the hits information
///             Date : oct/2016
///             Author : I. G. Irastorza
///
///_______________________________________________________________________________


#ifndef RestCore_TRestG4toHitsGetPhotoElectricProcess
#define RestCore_TRestG4toHitsGetPhotoElectricProcess

#include <TRestGas.h>
#include <TRestG4Event.h>
#include <TRestG4Metadata.h>
#include <TRestHitsEvent.h>

#include "TRestEventProcess.h"

class TRestG4toHitsGetPhotoElectricProcess:public TRestEventProcess {
    private:
#ifndef __CINT__
        TRestG4Event *fG4Event;//!
        TRestG4Metadata *fG4Metadata;//!
        TRestHitsEvent *fHitsEvent;//!

        vector <Int_t> fVolumeId;//!
#endif

        vector <TString> fVolumeSelection;

        void InitFromConfigFile();

        void Initialize();

        void LoadDefaultConfig();

    protected:

        //add here the members of your event process


    public:
        void InitProcess();
        void BeginOfEventProcess(); 
        TRestEvent *ProcessEvent( TRestEvent *eventInput );
        void EndOfEventProcess(); 
        void EndProcess();

        void LoadConfig( std::string cfgFilename, std::string name = "" );

        void PrintMetadata() { 

            BeginPrintProcess();

            for ( unsigned int n = 0; n < fVolumeSelection.size(); n++ )
                std::cout << "Volume added : " << fVolumeSelection[n] << std::endl;

            EndPrintProcess();

        }

        TString GetProcessName() { return (TString) "g4toHitsEventGetPhotoElectric"; }

        //Constructor
        TRestG4toHitsGetPhotoElectricProcess();
        TRestG4toHitsGetPhotoElectricProcess( char *cfgFileName );
        //Destructor
        ~TRestG4toHitsGetPhotoElectricProcess();

        ClassDef(TRestG4toHitsGetPhotoElectricProcess, 1);      // Transform a TRestG4Event event to a TRestHitsEvent (hits-collection event)
};
#endif

