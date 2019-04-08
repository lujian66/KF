//////////////////////////////////////////////////////////////////////////
///
///             RESTSoft : Software for Rare Event Searches with TPCs
///
///             TRestHitsKalmanFilterProcess.hpa
///
///              MAY 2019 : lujian
///
//////////////////////////////////////////////////////////////////////////


#ifndef RestCore_TRestHitsKalmanFilterProcess
#define RestCore_TRestHitsKalmanFilterProcess

#include <TRestHitsEvent.h>
#include "TRestEventProcess.h"
#include <Math/SVector.h>
typedef ROOT::Math::SVector<double,6>  SVector6;

class TRestHitsKalmanFilterProcess:public TRestEventProcess {

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

        double GetElectronVelocity(double E);
        double GetVelocityOffsetTheta02(double dZ, SVector6 XThcoli, Double_t E);

        void InitKalman();
        void Predictkalman();
        void Updatakalman();

        TString GetProcessName() { return (TString) "HitsKalmanFilter"; }

        //Constructor
        TRestHitsKalmanFilterProcess();
        TRestHitsKalmanFilterProcess( char *cfgFileName );
        //Destructor
        ~TRestHitsKalmanFilterProcess();

        ClassDef(TRestHitsKalmanFilterProcess, 1);      // Template for a REST "event process" class inherited from TRestEventProcess
};
#endif

