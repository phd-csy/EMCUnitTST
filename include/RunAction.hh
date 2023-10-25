#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class G4Timer;
/// Run action class

class RunAction : public G4UserRunAction
{
public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction(const G4Run *run) override;
    void EndOfRunAction(const G4Run *run) override;

private:
    G4Timer* timer;
};

#endif
