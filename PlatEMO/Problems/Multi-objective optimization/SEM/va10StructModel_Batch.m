classdef va10StructModel_Batch < PROBLEM
%<binary>
    % va10StructModel 结构方程模型多目标优化（工程改进版 · 方案4批处理）
    % 相对 va10StructModel 的改动（仅工程层面，目标函数定义不变）：
    %   1) 路径自包含：以本 .m 文件位置定位项目根目录，不再依赖 MATLAB 启动目录(pwd)
    %   2) 运行目录隔离：中间文件写入 ou\<问题名>\run_<进程ID>_<时间戳>\ ，并行互不覆盖
    %   3) 数据文件绝对路径写入 R 脚本，消除 R 工作目录依赖
    %   4) system 调用对路径加双引号，兼容含空格的目录名
    %   5) 方案4：每一代只调用一次 R —— 整代 N 个模型写进一个脚本，
    %      R 内 tryCatch 逐个拟合（单个体失败不影响整代），结果汇总写 res_gen<k>.csv；
    %      MATLAB 一次读回 N×8 矩阵，行号即个体编号
    %   6) 每代的 8 个原始指标 CSV 全部留档（res_gen<k>.csv），供日后重打分
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
                obj.D = 100;
            end
            obj.encoding = ones(1,obj.D)+3;
            obj.parameter.name = "va10StructModel_Batch";
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
            thisFile    = mfilename('fullpath');                      % ...\MMMOP\va10StructModel_Batch
            mmmopDir    = fileparts(thisFile);                        % ...\MMMOP
            projectRoot = fileparts(fileparts(fileparts(mmmopDir)));  % ...\<项目根，如 CopyModel-5>
            % ---- 运行目录隔离：进程 ID + 启动时间戳，进程内稳定（persistent）----
            persistent RUN_ID
            if isempty(RUN_ID)
                RUN_ID = sprintf('%d_%s', feature('getpid'), datestr(now,'yyyymmdd_HHMMSS'));
            end
            ouDir = fullfile(projectRoot,'ou',char(obj.parameter.name),['run_',RUN_ID]);
            if ~exist(ouDir,'dir')
                mkdir(ouDir);   % mkdir 会自动创建中间层目录
            end
            % ---- 代数计数：每代一批中间文件，全部留档 ----
            persistent GEN_CNT
            if isempty(GEN_CNT)
                GEN_CNT = 0;
            end
            GEN_CNT = GEN_CNT + 1;
            % ---- 数据文件使用绝对路径（R 中用正斜杠），消除 R 工作目录依赖 ----
            dataFileR = strrep(fullfile(mmmopDir,'va10.xlsx'),'\','/');

            %不可行解的大数惩罚（经 1-x / 3-x 变换后仍为极大正值，必被淘汰）
            PEN = 1e6;

            %用一个数组记录互惠关系
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
            disp("va10StructModel_Batch");

            %使用数据文件 va10.xlsx 中的数据
            Model_build1="model <- 'va1 =~ WR01A + WR01B + WR01C + WR01D";
            Model_build2="va2 =~ PEB02A + PEB02B + PEB02C + PEB02D + PEB02E";
            Model_build3="va3 =~ SI03A + SI03B + SI03C";
            Model_build4="va4 =~ CCB04A + CCB04B + CCB04C + CCB04D";
            Model_build5="va5 =~  EC05A + EC05B + EC05C + EC05D";
            Model_build6="va6 =~  EAT06A + EAT06B + EAT06C + EAT06D";
            Model_build7="va7 =~  PN07A + PN07B + PN07C";
            Model_build8="va8 =~ GTA08A + GTA08B + GTA08C + GTA08D";
            Model_build9="va9 =~ EK09A + EK09B + EK09C + EK09D";
            Model_build10="va10 =~ GPI10A + GPI10B + GPI10C";

            %单个体拟合块：tryCatch 容错，失败输出 NA 行，不中断整代
            fitBlock1='fit <- tryCatch(sem(model, data=datas), error=function(e) NULL)';
            fitBlock2='if (is.null(fit)) { a <- rep(NA,8) } else {';
            fitBlock3='a <- tryCatch(fitMeasures(fit,c("agfi","cfi","ifi","nnfi","nfi","rmsea","srmr","ecvi")), error=function(e) rep(NA,8)) }';
            fitBlock4='res <- rbind(res, a)';

            % ========== 生成整代的 R 脚本（方案4：一代一个脚本） ==========
            testFile = fullfile(ouDir,sprintf('test_gen%d.txt',GEN_CNT));
            outFile  = fullfile(ouDir,sprintf('out_gen%d.txt',GEN_CNT));   % R 控制台回显，调试用
            csvFile  = fullfile(ouDir,sprintf('res_gen%d.csv',GEN_CNT));   % N×8 拟合指标，留档
            csvFileR = strrep(csvFile,'\','/');
            %防止 R 失败时读到历史遗留的同名 CSV
            if exist(csvFile,'file')
                delete(csvFile);
            end
            fid=fopen(testFile,'wt');
            fprintf(fid,'%s\n', pre);
            fprintf(fid,'%s\n', pre1);
            fprintf(fid,'%s\n', pre2);
            fprintf(fid,'%s\n', 'res <- NULL');

            for each=1:N
                %测量模型
                fprintf(fid,'%s\n',Model_build1);
                fprintf(fid,'%s\n',Model_build2);
                fprintf(fid,'%s\n',Model_build3);
                fprintf(fid,'%s\n',Model_build4);
                fprintf(fid,'%s\n',Model_build5);
                fprintf(fid,'%s\n',Model_build6);
                fprintf(fid,'%s\n',Model_build7);
                fprintf(fid,'%s\n',Model_build8);
                fprintf(fid,'%s\n',Model_build9);
                fprintf(fid,'%s\n',Model_build10);

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

                %结束该个体的模型字符串，接拟合块
                fprintf(fid,'%s\n',"'");
                fprintf(fid,'%s\n', fitBlock1);
                fprintf(fid,'%s\n', fitBlock2);
                fprintf(fid,'%s\n', fitBlock3);
                fprintf(fid,'%s\n', fitBlock4);
            end  %结束N 结束each（仅生成脚本，不调R）

            %整代汇总写出
            fprintf(fid,'%s\n', 'colnames(res) <- c("agfi","cfi","ifi","nnfi","nfi","rmsea","srmr","ecvi")');
            fprintf(fid,'%s\n', sprintf('write.csv(res,"%s",row.names=FALSE)',csvFileR));
            fclose(fid);

            % ========== 整代只调用一次 R ==========
            %路径加双引号兼容含空格的目录名
            [status,~] = system(sprintf('R CMD BATCH --no-restore --no-save "%s" "%s"',testFile,outFile));

            %一次读回 N×8，行号即个体编号
            A = [];
            if status==0 && exist(csvFile,'file')
                A = readmatrix(csvFile);   % 跳过表头，读入 N×8 数值
            end

            for each=1:N
                objval = [];
                if ~isempty(A) && size(A,1)>=each && size(A,2)>=8 && ~any(isnan(A(each,1:8)))
                    a = A(each,1:8);
                    %目标组合与原定义一致：f1=agfi, f2=nnfi, f3=cfi+ifi+nfi+rmsea, f4=srmr, f5=ecvi
                    objval = [a(1), a(4), a(2)+a(3)+a(5)+a(6), a(7), a(8)];
                end
                if isempty(objval)
                    %不可行解大数惩罚：1-3目标经变换后仍为极大正值，4-5直接给极大正值
                    if obj.M<3
                        PopObj(each,1:obj.M)=-PEN;
                    else
                        PopObj(each,1:3)=-PEN;
                        PopObj(each,4:obj.M)=PEN;
                    end
                else
                    PopObj(each,:)=objval;
                end
            end
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
