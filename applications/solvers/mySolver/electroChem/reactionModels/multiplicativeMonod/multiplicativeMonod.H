#ifndef multiplicativeMonod_H
#define multiplicativeMonod_H

#include "reactionManager.H"

namespace Foam
{

class multiplicativeMonod : public reactionModel
{
    enum formType { firstOrder, monod };

    labelList substrateIndices_;
    List<formType> forms_;
    scalarList K_;
    scalarList cThreshold_;
    scalar value_; // qMax or k1, selected and dimensionally validated at startup
    scalar X_;
    scalar Xfilm_;
    const scalarField* filmIndicator_;

public:
    TypeName("multiplicativeMonod");

    multiplicativeMonod
    (
        const word& name,
        const dictionary& dict,
        const wordList& species
    );

    categoryType category() const override { return kinetic; }
    void setBiofilmIndicator(const scalarField* indicator) override
    { filmIndicator_ = indicator; }

    scalar kineticRate
    (
        const label celli,
        const scalarField& concentrations,
        const scalarField& nu
    ) const override;
};

} // End namespace Foam

#endif

