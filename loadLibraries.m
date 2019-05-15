% -- load missing libraries
function loadLibraries

loaded = javaclasspath('-dynamic');
cmFlag = 0;
orekitFlag = 0;
for k=1:length(loaded)
    if ~isempty(strfind(loaded{k}, 'math'))
        cmFlag = 1;
    end
    if ~isempty(strfind(loaded{k}, 'orekit'))
        orekitFlag = 1;
    end
end


if cmFlag == 0
    javaaddpath(fullfile(cd, 'lib', 'commons-math3-3.4.1.jar'));
end
if orekitFlag == 0
    javaaddpath(fullfile(cd, 'lib', 'orekit-7.0.jar'));
    javaaddpath(fullfile(cd, 'lib', 'orekit_custom-0.0.1-SNAPSHOT.jar'));
end

end

