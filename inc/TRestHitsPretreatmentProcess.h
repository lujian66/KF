//////////////////////////////////////////////////////////////////////////
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestHitsPretreatmentProcess.hpa
///
///              MAY 2019 : lujian
///
//////////////////////////////////////////////////////////////////////////


#ifndef RestCore_TRestHitsPretreatmentProcess
#define RestCore_TRestHitsPretreatmentProcess

#include <TRestHitsEvent.h>
#include "TRestEventProcess.h"

class TRestHitsPretreatmentProcess:public TRestEventProcess {

    private:

#ifndef __CINT__
        TRestHitsEvent *fInputHitsEvent;//!
        TRestHitsEvent *fOutputHitsEvent;//!
#endif

        void InitFromConfigFile();

        void Initialize();

    protected:

     /*   Double_t fStartingDistance;
        Double_t fMinimumDistance;
        Double_t fDistanceFactor;
        Double_t fMaxNodes;
        Double_t ZMax;
        Double_t ZMin;
	*/

    public:
        void InitProcess();
        void BeginOfEventProcess(); 
        TRestEvent *ProcessEvent( TRestEvent *eventInput );
        void EndOfEventProcess(); 
        void EndProcess();
        void LoadDefaultConfig( );

        void LoadConfig( std::string cfgFilename, std::string name = "" );

        void PrintMetadata() { 

            BeginPrintProcess();

         /*   std::cout << " Starting distance : " << fStartingDistance << std::endl;
            std::cout << " Minimum distance : " << fMinimumDistance << std::endl;
            std::cout << " Distance step factor : " << fDistanceFactor << std::endl;
            std::cout << " Maximum number of nodes : " << fMaxNodes << std::endl;
*/
            EndPrintProcess();
        }

        TString GetProcessName() { return (TString) "HitsPretreatment"; }

        //Constructor
        TRestHitsPretreatmentProcess();
        TRestHitsPretreatmentProcess( char *cfgFileName );
        //Destructor
        ~TRestHitsPretreatmentProcess();

        ClassDef(TRestHitsPretreatmentProcess, 1);      // Template for a REST "event process" class inherited from TRestEventProcess
};
#endif

