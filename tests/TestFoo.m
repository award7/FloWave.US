classdef TestFoo < matlab.unittest.TestCase
    methods(Test)
        function test_foo(testCase)
            expected = "Hello, world";
            actual = foo();
            testCase.verifyEqual(actual, expected);
        end
    end
end
