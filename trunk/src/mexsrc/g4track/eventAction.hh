#include "G4UserEventAction.hh"
#include "globals.hh"

/// Event action class

class eventAction : public G4UserEventAction
{
  public:
    eventAction();
    virtual ~eventAction();

    virtual void  BeginOfEventAction(const G4Event* );
    virtual void    EndOfEventAction(const G4Event* );
};


