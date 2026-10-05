classdef va5StructModel_New < PROBLEM
%<binary>
    % va5 结构方程模型多目标优化（工程改进版）
    % 相对 va3StructModel 的改动（仅工程层面，目标函数定义不变）：
    %   1) 路径自包含：以本 .m 文件位置定位项目根目录，不再依赖 MATLAB 启动目录(pwd)
    %   2) 运行目录隔离：中间文件写入 ou\<问题名>\run_<进程ID>_<时间戳>\ ，
    %      不同算法 / 不同 MATLAB 实例并行运行同一问题时互不覆盖，且历史运行可追溯
    %   3) 数据文件绝对路径写入 R 脚本，消除 R 工作目录依赖
    %   4) system 调用对路径加双引号，兼容含空格的目录名
    %   5) R 端 8 个拟合指标写 CSV，MATLAB 端 readmatrix 直读，
    %      废掉按"倒数第 5 行"解析的脆弱逻辑；CSV 留档可供日后重打分
    %   6) 检查 system 返回值与 CSV 是否生成，R 失败立即判不可行
    %   7) 不可行惩罚由 inf 改为大数 PEN（1e6），避免 HV/IGD/归一化出 NaN

    properties(Access = private)
    end
    methods
        %% Initialization
        function Setting(obj)
            % Parameter setting
            if isempty(obj.M)
                obj.M = 5;
            end
            if isempty(obj.D)
                obj.D = 25;
            end
            obj.encoding = ones(1,obj.D)+3;
            obj.parameter.name = "va5StructModel_New";
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
        end
        %% Repair infeasible solutions
        function PopDec = CalDec(obj,PopDec)
            N = size(PopDec,1);
            VarNum = sqrt(size(PopDec,2));
            for t = 1:N
                for p = 1: VarNum
                    PopDec(t,VarNum*(p-1)+p) = 0;
                end
            end
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            % ---- 路径解析：以本文件位置为准，与 MATLAB 当前目录无关 ----
            % 本文件位于 <项目根>\Problems\Multi-objective optimization\MMMOP\
            thisFile    = mfilename('fullpath');                    % ...\MMMOP\va5StructModel_New
            mmmopDir    = fileparts(thisFile);                      % ...\MMMOP
            projectRoot = fileparts(fileparts(fileparts(mmmopDir)));% ...\<项目根，如 CopyModel-5>
            % ---- 运行目录隔离：进程 ID + 启动时间戳，进程内稳定（persistent）----
            persistent RUN_ID
            if isempty(RUN_ID)
                RUN_ID = sprintf('%d_%s', feature('getpid'), datestr(now,'yyyymmdd_HHMMSS'));
            end
            ouDir = fullfile(projectRoot,'ou',char(obj.parameter.name),['run_',RUN_ID]);
            if ~exist(ouDir,'dir')
                mkdir(ouDir);   % mkdir 会自动创建中间层目录
            end
            % ---- 数据文件使用绝对路径（R 中用正斜杠），消除 R 工作目录依赖 ----
            dataFileR = strrep(fullfile(mmmopDir,'va5.xlsx'),'\','/');

            %不可行解的大数惩罚（经 1-x / 3-x 变换后仍为极大正值，必被淘汰）
            PEN = 1e6;

            %用一个数组记录互惠关系
            %在这个函数里面计算目标函数值，输入是一个N*D的矩阵每一行是一个个体的决策向量，要输出一个目标函数值矩阵N*M 分别和
            %决策向量每一行对应
            N=size(PopDec,1);
            VarNum=sqrt(size(PopDec,2));%潜变量数目
            relationMap=cell(N,1);
            for Flag_Pop=1:N
                relation=[];
                %第i个潜变量找到第j个潜变量如果是互惠关系的话就先设置为没有关系然后保存互惠关系的矩阵
                for i=1:VarNum
                    for j=1:VarNum
                        if and(PopDec(Flag_Pop,VarNum*(i-1)+j),PopDec(Flag_Pop,VarNum*(j-1)+i))
                            PopDec(Flag_Pop,VarNum*(i-1)+j)=0;
                            PopDec(Flag_Pop,VarNum*(j-1)+i)=0;
                            relation=[relation;i,j];
                        end
                    end
                end
                relationMap{Flag_Pop,1}=relation;
            end
            %建立目标值矩阵固定5个目标
            PopObj=ones(N,5);
            %对结构模型进行建模
            pre="library(lavaan)";
            pre1="library(openxlsx)";
            pre2=sprintf('datas<-read.xlsx("%s")',dataFileR);
            disp("va5_New");

            %使用va5.xlsx中的数据
            Model_build1="model <- 'va1 =~ q1 + q2+q3+q4 ";
            Model_build2="va2 =~ q5 + q6+q7+q8+q9 ";
            Model_build3="va3 =~  q10 + q11+q12+q13+q14 ";
            Model_build4="va4 =~ q15 + q16 + q17+q18+q19+q20+q21+q22+q23 ";
            Model_build5="va5 =~ q24  ";

            %agfi 和 gfi 的上限都是1 而且是越大越好。rmr要求<0.05同样越小越好。但是是不是0.05的时候就可以当做可忽略的对象呢？
            cac="fit <- sem(model, data =datas)";
            cac_valid= "a<-fitMeasures(fit,c(""agfi"",""cfi"",""ifi"",""nnfi"",""nfi"",""rmsea"",""srmr"",""ecvi""))";

            %对每一个个体进行分析流程
            for each=1:N
                testFile = fullfile(ouDir,'test.txt');
                outFile  = fullfile(ouDir,sprintf('ou%d.txt',each));   % R 控制台回显，调试用
                csvFile  = fullfile(ouDir,sprintf('res%d.csv',each));  % 8 个拟合指标
                csvFileR = strrep(csvFile,'\','/');
                %防止 R 失败时读到上一代残留的同名 CSV
                if exist(csvFile,'file')
                    delete(csvFile);
                end
                fid=fopen(testFile,'wt');
                fprintf(fid,'%s\n', pre);

                %使用外来数据
                fprintf(fid,'%s\n', pre1);
                fprintf(fid,'%s\n', pre2);
                fprintf(fid,'%s\n',Model_build1);
                fprintf(fid,'%s\n', Model_build2);
                fprintf(fid,'%s\n',Model_build3);
                fprintf(fid,'%s\n',Model_build4);
                fprintf(fid,'%s\n',Model_build5);

                %对每一个个体进行生成结构模型代码
                for  i=1:VarNum
                    %对each个体的第i个潜变量建模
                    temp=char(sprintf("va%d~",i));
                    flag=0;
                    for j=1:VarNum
                        %判断其他个体对第i个潜变量的影响
                        if PopDec(each,VarNum*(i-1)+j)==1
                            if flag==0
                                temp2=char(sprintf("va%d",j));
                                temp=[temp,temp2];
                                flag=flag+1;
                            else
                                temp2=char(sprintf("+va%d",j));
                                temp=[temp,temp2];
                            end
                        end
                    end
                    %如果有影响就输入文件
                    if and(length(temp)>4,i<10)
                        fprintf(fid,'%s\n',string(temp));
                    elseif and(length(temp)>5,i>=10)
                        fprintf(fid,'%s\n',string(temp));
                    end
                end

                %考虑互惠关系
                relationCell=relationMap{each,1};
                if size(relationCell,1)>=1
                    for reci_size=1:size(relationCell,1)
                        Array=relationCell(reci_size,:);
                        PopDec(each,VarNum*(Array(1)-1)+Array(2))=1;
                        PopDec(each,VarNum*(Array(2)-1)+Array(1))=1;
                        temp=char(sprintf("va%d~~va%d",Array(1),Array(2)));
                        fprintf(fid,'%s\n',string(temp));
                    end
                end

                %使用自己的数据
                fprintf(fid,'%s\n',"'");
                fprintf(fid,'%s\n', cac);
                fprintf(fid,'%s\n',cac_valid);
                %8 个原始拟合指标写 CSV（留档，供日后重打分 / 换目标组合）
                fprintf(fid,'%s\n', sprintf('write.csv(t(a),"%s",row.names=FALSE)',csvFileR));
                fclose(fid);
                %调用R语言构建目标函数的解，路径加双引号兼容含空格的目录名
                [status,~] = system(sprintf('R CMD BATCH --no-restore --no-save "%s" "%s"',testFile,outFile));
                %status 非 0（R 未启动/崩溃）或 CSV 未生成（lavaan 不收敛/报错）均判不可行
                objval = [];
                if status==0 && exist(csvFile,'file')
                    a = readmatrix(csvFile);   % 跳过表头，读入 1x8 数值行
                    if numel(a)>=8 && ~any(isnan(a(1:8)))
                        %目标组合与原定义一致：f1=agfi, f2=nnfi, f3=cfi+ifi+nfi+rmsea, f4=srmr, f5=ecvi
                        objval = [a(1), a(4), a(2)+a(3)+a(5)+a(6), a(7), a(8)];
                    end
                end
                %目前想法是把那些不符合要求的结果进行最大化处理最后会被淘汰掉
                %并且把要做相同或者类似处理的目标值放在相邻的地方
                if isempty(objval)
                    %如果是不合适的解就对其进行一些处理，1-3是越大越好要求最小值
                    if obj.M<3
                        PopObj(each,1:obj.M)=-PEN;
                    else
                        %先把1：3变成负的，被一减了之后就变成了正的。越大越容易被淘汰
                        PopObj(each,1:3)=-PEN;
                        PopObj(each,4:obj.M)=PEN;
                    end
                else
                    PopObj(each,:)=objval;
                end
            end  %结束N 结束each
            for i=1:obj.M
                if i<=2
                    PopObj(:,i)=1-PopObj(:,i);
                elseif i==3
                    PopObj(:,i)=3-PopObj(:,i);
                end
            end
            disp(PopObj);
        end
    end
end
