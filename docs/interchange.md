# Interchange and ownership

Input and output use explicit schema discriminants. Generated schemas live in
`schemas/`; the Rust model is authoritative. Unknown semantic fields, recursive
duplicate keys, positional substitutes for objects, unsupported schemas, and
non-finite values are rejected. New semantic fields require negotiated schema
revision. Object ordering and observation ordering are canonicalized by opaque
string identities; bulk endpoint integers are transport indexes only.

Input requires frame, mode, stations, observations, controls, relationships, and
algorithm. Frame provenance and observation reading identities default to empty.
Spatial mode has no preserved-Z map; horizontal mode requires a complete map.
Hard/soft controls supply one value per listed active axis. Soft covariance is
square in that axis order. Recorded relationships retain evidence without
creating a numerical constraint. Finite ties are ordinary displacement factors.

`tests/fixtures/controlled_tree.json` has two measured ten-metre legs between
controls twenty-one metres apart. Standard deviations one and two yield B at
10.2 metres, corrections 0.2 and 0.8, objective 0.2, no graph cycle, and scalar
redundancy three. Expected values are analytic, independent of the implementation.

Preparation owns copied validated state. Source mutation after preparation does
not change a solve. Source mutation during buffer copying is prohibited; copying
is not an atomic snapshot. Read-only memoryviews can still have writable aliases.
A result must only be applied to its associated preparation and source baseline.
