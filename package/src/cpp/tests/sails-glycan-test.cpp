//
// Created by Jordan Dialpuri on 31/07/2024.
//

#include "../../include/sails-glycan.h"

#include <gtest/gtest.h>

TEST(GlycanTest, TestFunctionality) {
    Sails::Glycan glycan = {};
    Sails::Glycosite g1 = {0,0,0};
    Sails::Glycosite g2 = {0,0,1};

    EXPECT_TRUE(glycan.empty()); // glycan should start empty

    EXPECT_EQ(glycan.begin(), glycan.end()); // iterators should be the same with an empty glycan

    glycan.add_sugar("A", 1, g1); // add first sugar
    EXPECT_EQ(glycan.size(), 1);
    EXPECT_NE(glycan.sugars[g1], nullptr);

    glycan.add_sugar("A", 1, g1); // add first sugar again
    EXPECT_EQ(glycan.size(), 1);

    glycan.add_sugar("A", 1, g2); // add second sugar
    EXPECT_EQ(glycan.size(), 2);
    EXPECT_NE(glycan.sugars[g2], nullptr);

    glycan.add_linkage(g1, g2, "A", "B"); // add linkage
    EXPECT_EQ(glycan.adjacency_list.size(), 1);

    auto s1 = glycan.sugars[g1].get();
    auto s2 = *glycan.adjacency_list[s1].begin();
    EXPECT_NE(s1, nullptr);
    EXPECT_NE(s2, nullptr);

    auto prev = glycan.find_previous_sugar(s2);
    EXPECT_TRUE(prev.has_value());
    EXPECT_EQ(s1, prev.value());

    auto linked = glycan.remove_sugar(s2, false);
    EXPECT_EQ(s1, linked);
}

TEST(SugarTest, TestFunctionality) {
    Sails::Glycosite g1 = {0,0,0};
    Sails::Glycosite g2 = {0,0,1};

    Sails::Sugar s1 = {"A", 0, g1};
    Sails::Sugar s2 = {"B", 1, g2};

    EXPECT_LT(s1, s2);
}
