classdef NonRegTestClass < matlab.unittest.TestCase
    properties
    end
    methods (Test,TestTags = {'All'})
        function Init(testCase)
            a = which(mfilename);
            a = regexp(a,'@','split');
            a = a{1};
            addpath(a);
        end
    end
    methods (Test,TestTags = {'All','Elastic'})
        function ElasticAnisotropic2D(testCase)
            AnisotropicElastic2D
        end
        function ElasticAnisotropic3D(testCase)
            AnisotropicElastic3D
        end
        function ElasticIsotropic2D(testCase)
            IsotropicElastic2D
        end
        function ElasticIsotropic3D(testCase)
            IsotropicElastic3D
        end
    end
    methods (Test,TestTags = {'All','Acoustic'})
        function AcousticAnisotropic2D(testCase)
            AnisotropicAcoustic2D
        end
        function AcousticAnisotropic3D(testCase)
            AnisotropicAcoustic3D
        end
        function AcousticIsotropic2D(testCase)
            IsotropicAcoustic2D
        end
        function AcousticIsotropic3D(testCase)
            IsotropicAcoustic3D
        end
        function AcousticIsotropic3DDifCordSystems(testCase)
            IsotropicAcoustic3DDifCordSystems
        end
        function AcousticAnisotropic3DDifCordSystems(testCase)
            AnisotropicAcoustic3DDifCordSystems
        end
    end
    methods (Test,TestTags = {'2D'})
        function Acoustic2DIsotropic(testCase)
            IsotropicAcoustic2D
        end
        function Acoustic2DAnisotropic(testCase)
            AnisotropicAcoustic2D
        end
        function Elastic2DAnisotropic(testCase)
            AnisotropicElastic2D
        end
        function Elastic2DIsotropic(testCase)
            IsotropicElastic2D
        end
    end
    methods (Test,TestTags = {'3D'})
        function Acoustic3DIsotropic(testCase)
            IsotropicAcoustic3D
        end
        function Acoustic3DAnisotropic(testCase)
            AnisotropicAcoustic3D
        end
        function Elastic3DAnisotropic(testCase)
            AnisotropicElastic3D
        end
        function Elastic3DIsotropic(testCase)
            IsotropicElastic3D
        end
        function Acoustic3DDifCordSystemsIsotropic(testCase)
            IsotropicAcoustic3DDifCordSystems
        end
        function Acoustic3DDifCordSystemsAnisotropic(testCase)
            AnisotropicAcoustic3DDifCordSystems
        end
    end
    methods (Test,TestTags = {'All','Material'})

    end
end