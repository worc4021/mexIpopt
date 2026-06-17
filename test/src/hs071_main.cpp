// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10

// [MAIN]
#include "IpIpoptApplication.hpp"
#include "hs071_nlp.hpp"
#include <gtest/gtest.h>
#include <iostream>
#include <map>


class IpoptHS071 
   : public ::testing::Test 
{
   public:
      std::map<std::string,double> realOptions;
      std::map<std::string,int> intOptions;
      std::map<std::string,std::string> strOptions;
   private:
      Ipopt::SmartPtr<Ipopt::TNLP> mynlp;
      Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

   public:
      IpoptHS071()
         : mynlp(new HS071_NLP()),
           app(IpoptApplicationFactory())
      {
      }

      void TestBody() override {
         for (const auto& [key, value] : realOptions) {
            EXPECT_TRUE(setOption(key, value))
               << "Failed to set real option '" << key << "'"  ;
         }
         for (const auto& [key, value] : intOptions) {
            EXPECT_TRUE(setOption(key, value))
               << "Failed to set integer option '" << key << "'"  ;
         }
         for (const auto& [key, value] : strOptions) {
            EXPECT_TRUE(setOption(key, value))
               << "Failed to set string option '" << key << "'"  ;
         }

         ApplicationReturnStatus status;
         status = app->Initialize();
         EXPECT_EQ( status, Solve_Succeeded ) << "Error during initialization!";

         // Ask Ipopt to solve the problem
         status = app->OptimizeTNLP(mynlp);

         EXPECT_EQ( status, Solve_Succeeded ) << "The problem FAILED!";
      }

      bool setOption(const std::string& key, double value) {
         return app->Options()->SetNumericValue(key, value);
      }
      bool setOption(const std::string& key, int value) {
         return app->Options()->SetIntegerValue(key, value);
      }
      bool setOption(const std::string& key, const std::string& value) {
         return app->Options()->SetStringValue(key, value);
      }
};

TEST_F(IpoptHS071, MA27)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "ma27";
   test.TestBody();
}

TEST_F(IpoptHS071, MA57)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "ma57";
   test.TestBody();
}

TEST_F(IpoptHS071, MA77)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "ma77";
   test.TestBody();
}

TEST_F(IpoptHS071, MA86)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "ma86";
   test.TestBody();
}

TEST_F(IpoptHS071, MA97)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "ma97";
   test.TestBody();
}

TEST_F(IpoptHS071, MUMPS)
{
   IpoptHS071 test;
   test.strOptions["linear_solver"] = "mumps";
   test.TestBody();
}


int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}