#include "definitions.h"

Atom::Atom(Vector v, Metal metal) {
    coord = v;
    type = metal;
}

bool Atom::operator==(Atom other) const {
    return this->coord == other.coord && this->type == other.type;
}