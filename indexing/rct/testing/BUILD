package(
    default_visibility = ["//visibility:public"],
)

cc_test(
    name = "dense_feature_vector_test",
    srcs = ["dense_feature_vector_test.cc"],
    size = "small",
    deps = [
        "//indexing/rct/public:dense_feature_vector",
        "//third_party/gtest:gtest_main",
    ],
)

cc_test(
    name = "metrics_test",
    srcs = ["metrics_test.cc"],
    size = "small",
    deps = [
        "//indexing/rct/public:dense_feature_vector",
        "//indexing/rct/public:metrics",
        "//third_party/gtest:gtest_main",
    ],
)

cc_test(
    name = "rank_cover_tree_test",
    srcs = ["rank_cover_tree_test.cc"],
    size = "small",
    deps = [
        "//indexing/rct/public:dense_feature_vector",
        "//indexing/rct/public:metrics",
        "//indexing/rct/public:rank_cover_tree",
        "//third_party/gtest:gtest_main",
    ],
)
