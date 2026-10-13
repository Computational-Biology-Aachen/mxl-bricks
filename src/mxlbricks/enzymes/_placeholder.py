# Copyright 2025 Marvin van Aalst, Anna Matuszyńska

# x1 = default_name(x1, n.)
# x2 = default_name(x2, n.)
# x3 = default_name(x3, n.)
# x4 = default_name(x4, n.)
# x5 = default_name(x5, n.)
# x6 = default_name(x6, n.)


# def add_rxn_1_1(
#     *,
#     rxn: str | None = None,
#     s1: str | None = None,
#     p1: str | None = None,
# ) -> Model:
#     rxn = default_name(rxn, n.)
#     s1 = default_name(s1, n.)
#     p1 = default_name(p1, n.)
#     model.add_reaction(
#         name=rxn,
#         stoichiometry=filter_stoichiometry(
#             model,
#             {
#                 s1: -1,
#                 p1: 1,
#             },
#         ),
#         args=[
#             s1,
#             p1,
#             default_vmax(
#                 model,
#                 rxn=rxn,
#                 e0=e0,
#                 kcat=kcat,
#                 e0_default=1.0,  # Source
#                 kcat_default=None,  # Source
#             ),
#             default_kms(model, rxn=rxn, name=kms, default=None),
#             default_kmp(model, rxn=rxn, name=kmp, default=None),
#             default_keq(model, rxn=rxn, name=keq, default=None),
#         ],
#     )

#     return model

# def add_rxn_2_2(
#     *,
#     rxn: str | None = None,
#     x1: str | None = None,
#     x2: str | None = None,
#     x3: str | None = None,
#     x4: str | None = None,
# ) -> Model:
#     rxn = default_name(rxn, n.)
#     x1 = default_name(x1, n.)
#     x2 = default_name(x2, n.)
#     x3 = default_name(x3, n.)
#     x4 = default_name(x4, n.)
#     model.add_reaction(
#         name=rxn,
#         stoichiometry=filter_stoichiometry(
#             model,
#             {
#                 x1: -1,
#                 x2: -1,
#                 x3: 1,
#                 x4: 1,
#             },
#         ),
#         args=[
#             x1,
#             x2,
#             x3,
#             x4,
#             default_vmax(
#                 model,
#                 rxn=rxn,
#                 e0=e0,
#                 kcat=kcat,
#                 e0_default=1.0,  # Source
#                 kcat_default=None,  # Source
#             ),
#             default_kms(model, rxn=rxn, name=kms, default=None),
#             default_kmp(model, rxn=rxn, name=kmp, default=None),
#             default_keq(model, rxn=rxn, name=keq, default=None),
#         ],
#     )

#     return model


# def add_rxn_3_3(
#     *,
#     rxn: str | None = None,
#     s1: str | None = None,
#     s2: str | None = None,
#     s3: str | None = None,
#     p1: str | None = None,
#     p2: str | None = None,
#     p3: str | None = None,
# ) -> Model:
#     rxn = default_name(rxn, n.)
#     s1 = default_name(s1, n.)
#     s2 = default_name(s2, n.)
#     s3 = default_name(s3, n.)
#     p1 = default_name(p1, n.)
#     p2 = default_name(p2, n.)
#     p3 = default_name(p3, n.)
#     model.add_reaction(
#         name=rxn,
#         stoichiometry=filter_stoichiometry(
#             model,
#             {
#                 s1: -1,
#                 s2: -1,
#                 s3: -1,
#                 p1: 1,
#                 p2: 1,
#                 p3: 1,
#             },
#         ),
#         args=[
#             s1,
#             s2,
#             s3,
#             p1,
#             p2,
#             p3,
#             default_vmax(
#                 model,
#                 rxn=rxn,
#                 e0=e0,
#                 kcat=kcat,
#                 e0_default=1.0,  # Source
#                 kcat_default=None,  # Source
#             ),
#             default_kms(model, rxn=rxn, name=kms, default=None),
#             default_kmp(model, rxn=rxn, name=kmp, default=None),
#             default_keq(model, rxn=rxn, name=keq, default=None),
#         ],
#     )

#     return model


# def add_kre_rxn_2_1(
#     *,
#     rxn: str | None = None,
#     s1: str | None = None,
#     s2: str | None = None,
#     p1: str | None = None,
# ) -> Model:
#     rxn = default_name(rxn, n.)
#     s1 = default_name(s1, n.)
#     s2 = default_name(s2, n.)
#     p1 = default_name(p1, n.)
#     model.add_reaction(
#         name=rxn,
#         stoichiometry=filter_stoichiometry(
#             model,
#             {
#                 s1: -1,
#                 s2: -1,
#                 s3: -1,
#                 p1: 1,
#                 p2: 1,
#                 p3: 1,
#             },
#         ),
#         args=[
#             s1,
#             s2,
#             p1,
#             default_kre(model, rxn=rxn, name=kre, default=800000000.0),
#             default_keq(model, rxn=rxn, name=keq, default=None),
#         ],
#     )

#     return model
