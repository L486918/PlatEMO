classdef va6StructModel < PROBLEM
%<binary>
    
    %------------------------------- Reference --------------------------------
    % E. Zitzler and L. Thiele, Multiobjective evolutionary algorithms: A
    % comparative case study and the strength Pareto approach, IEEE
    % Transactions on Evolutionary Computation, 1999, 3(4): 257-271.
    %------------------------------- Copyright --------------------------------
    % Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
    % research purposes. All publications which use this platform or any code
    % in the platform should acknowledge the use of "PlatEMO" and reference "Ye
    % Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
    % for evolutionary multi-objective optimization [educational forum], IEEE
    % Computational Intelligence Magazine, 2017, 12(4): 73-87".
    %--------------------------------------------------------------------------
    
    
    %这个问题可能会需要测量模型的方差%
    %在这个文件导入对应的数据并计算%
    
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
                obj.D = 36;
            end
            obj.encoding = ones(1,obj.D)+3;
            obj.parameter.name="va6StructModel";
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
            % 笔者路径为D:\Codes\StructuralEquationModel\算法名\Problems\Multi-objectivptimization\MMMOP
            % 提取D:\Codes\StructuralEquationModel\算法名\ou 作为中间文件存储路径
            currentPath=pwd;
            parts=strsplit(currentPath,'\');
            parentPath=strjoin(parts(1:4),'\');            
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
            %建立目标值矩阵固定8个目标
            PopObj=ones(N,5);
            %对结构模型进行建模
            pre="library(lavaan)";
            pre1="library(openxlsx)";
            pre2="datas<-read.xlsx(""va6.xlsx"")";
            disp("va6");
            
            
            
            Model_build1="model <- 'va1 =~ STM01A + STM01B + STM01C + STM01D";  %va1即为StM
            Model_build2="va2 =~ ASE02A + ASE02B + ASE02C + ASE02D + ASE02E";  %va2即为ASE
            Model_build3="va3 =~ TCS03A + TCS03B + TCS03C + TCS03D + TCS03E ";  %TcS
            Model_build4="va4 =~ PEI04A + PEI04B + PEI04C + PEI04D ";  %PeI
            Model_build5="va5 =~ LES05A + LES05B + LES05C + LES05D + LES05E + LES05F "; %LeS
            Model_build6="va6 =~ ACP06A + ACP06B + ACP06C + ACP06D + ACP06E";  %AcP
            
            
            %agfi 和 gfi 的上限都是1 而且是越大越好。rmr要求<0.05同样越小越好。但是是不是0.05的时候就可以当做可忽略的对象呢？
            cac="fit <- sem(model, data =datas)";
            cac_std="fit <- sem(model, data =datas,std.lv=T)";
            cac_valid= "a<-fitMeasures(fit,c(""agfi"",""cfi"",""ifi"",""nnfi"",""nfi"",""rmsea"",""srmr"",""ecvi""))";
            c= "i<-c(a[1],a[4],a[2]+a[3]+a[5]+a[6],a[7],a[8])";
            e="i";
            %如果目标数量超过可以评估的数目就使用cac1输出全部
            %结构模型需要进行分析后构建   StructModel1 StructModel2 StructModel3 变量的数目和潜变量的数目相等
            %潜变量用va表示
            %先把string变成char
            %可以通过向char添加新的串达到目的
            
            %对每一个个体进行分析流程
            %  for each=1:N
            %  输入测量模型
            % for  i=1:VarNum
            % 针对每一个潜变量i
            %for j=1:VarNum
            %判断所有潜变量j对i的影响并且输入
            %结束for j 潜变量i的构建已经完成
            %如果构建后的潜变量如果满足就输入文件
            %结束for i
            %运行代码
            %更新目标值矩阵
            %结束
            
            for each=1:N
 
                fid=fopen(sprintf('%s\\ou\\%s\\test.txt',parentPath,obj.parameter.name),'wt');
                fprintf(fid,'%s\n', pre);
              
                
                %使用外来数据
                fprintf(fid,'%s\n', pre1);
                fprintf(fid,'%s\n', pre2);
                fprintf(fid,'%s\n',Model_build1);
                fprintf(fid,'%s\n', Model_build2);
                fprintf(fid,'%s\n',Model_build3);
                fprintf(fid,'%s\n',Model_build4);
                fprintf(fid,'%s\n',Model_build5);
                fprintf(fid,'%s\n',Model_build6);
                
                
                
                %对每一个个体进行生成结构模型代码
                for  i=1:VarNum
                    %对each个体的第i个潜变量建模
                    temp=char(sprintf("va%d~",i)); % 果 i 
                    flag=0;
                    for j=1:VarNum
                        %判断其他个体对第i个潜变量的影响
                        if PopDec(each,VarNum*(i-1)+j)==1 % 因 j
                            if flag==0
                                temp2=char(sprintf("va%d",j));
                                temp=[temp,temp2];
                                flag=flag+1;
                            else
                                %说明这个潜变量不仅仅影响一个，还影响多个
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
                fprintf(fid,'%s\n',c);
                fprintf(fid,'%s\n',e);
                %fprintf(fid,'%s\n',sum);
                fclose(fid);
                %到目前为止要输入的模型已经构建完成接下来要做的事调用R语言然后构建目标函数的解然后读取out.txt文件选取要用到的指标
                system(sprintf('R CMD BATCH --no-restore --no-save %s\\ou\\%s\\test.txt %s\\ou\\%s\\ou%d.txt',parentPath,obj.parameter.name,parentPath,obj.parameter.name,each));
                lines= get_lines(sprintf('%s\\ou\\%s\\ou%d.txt',parentPath,obj.parameter.name,each));
                objval=get_objval(sprintf('%s\\ou\\%s\\ou%d.txt',parentPath,obj.parameter.name,each),lines);
                try
                   objval=str2num(objval);
                catch
                    objval=[];
                end
                %目前想法是把那些不符合要求的结果进行最大化处理最后会被淘汰掉
                %并且把要做相同或者类似处理的目标值放在相邻的地方
                if isempty(objval)
                    %如果是不合适的解就对其进行一些处理，1-3是越大越好要求最小值
                    if obj.M<3
                        PopObj(each,1:obj.M)=-inf;
                    else
                        %先把1：5变成负的，被一减了之后就变成了正的。越大越容易被淘汰
                        PopObj(each,1:3)=-inf;
                        PopObj(each,4:obj.M)=inf;
                    end
                else
                    PopObj(each,:)=objval;
                end
            end  %结束N 结束each
            if all(all(isinf(PopObj)))
                error("检测到有全部是inf的矩阵，检查是否存在问题")
            end
            for i=1:obj.M
                if i<=2
                    PopObj(:,i)=1-PopObj(:,i);
                elseif i==3
                    PopObj(:,i)=3-PopObj(:,i);
                end
            end
        end
    end
end