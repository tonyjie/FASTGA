# FastGA 命令选项代码位置说明

本文档指出 FastGA 各个命令选项在源代码中的定义位置。

## 文件位置
主要源代码文件：`FastGA.c`

---

## 1. 命令选项定义位置

### 选项声明和初始化

**位置：** `FastGA.c` 第 62-86 行

```62:86:FastGA.c
static char *Usage[] = { "[-vkMS] [-L:<log:path>] [-T<int(8)>] [-P<dir($TMPDIR)>] [<format(-paf)>]",
                         "[-f<int(10)>] [-c<int(85)> [-s<int(1000)>] [-l<int(100)>] [-i<float(.7)]",
                         "<source1:path>[<precursor>] [<source2:path>[<precursor>]]"
                       };

static int    NEW_GIX;     //  Is this a rev 2.0 style GIX?

static int    FREQ;        //  -f: Adaptamer frequency cutoff parameter
static int    VERBOSE;     //  -v: Verbose output
static char  *LOG_PATH;    //  -L: Log file path name
static FILE  *LOG_FILE;    //      Log file handle
static int    CHAIN_BREAK; //  -s
static int    CHAIN_MIN;   //  -c
static int    ALIGN_MIN;   //  -a
static double ALIGN_RATE;  //  1.-e
static int    NTHREADS;    //  -T
static char  *SORT_PATH;   //  -P
static int    KEEP;        //  -k
static int    SYMMETRIC;   //  -S
static int    SOFT_MASK;   //  -M
static int    SELF;        //  Comparing A to A, or A to B?
static int    OUT_TYPE;    //  -paf = 0; -psl = 1; -one = 2
static int    OUT_OPT;     //  -paf = 0, -m = PAFM; -x = PAFX; -s = PAFS; -S = PAFL
static char  *ONE_PATH;    //  -one option path
static char  *ONE_ROOT;    //  -one option path
```

---

## 2. 命令行参数解析

**位置：** `FastGA.c` 第 4433-4565 行（main函数中的参数解析部分）

### 2.1 默认值初始化

```4448:4463:FastGA.c
    FREQ = 10;
    CHAIN_BREAK = 2000;   //  2x in anti-diagonal space
    CHAIN_MIN   =  170;
    ALIGN_MIN   =  100;
    ALIGN_RATE  = .3;
    SORT_PATH = getenv("TMPDIR");
    if (SORT_PATH == NULL)
      SORT_PATH = ".";
    NTHREADS    = 8;

    OUT_TYPE    = 0;
    OUT_OPT     = 0;
    ONE_PATH    = NULL;
    ONE_ROOT    = NULL;
    LOG_PATH    = NULL;
    LOG_FILE    = NULL;
```

### 2.2 参数解析循环

```4466:4565:FastGA.c
    j = 1;
    for (i = 1; i < argc; i++)
      if (argv[i][0] == '-')
        switch (argv[i][1])
        { default:
            ARG_FLAGS("vkMS")
            break;
          case '1':
            if (strncmp(argv[i]+1,"1:",2) == 0)
              { OUT_TYPE = 2;
                ONE_PATH = PathTo(argv[i]+3);
                ONE_ROOT = Root(argv[i]+3,".1aln");
                test = fopen(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"),"w");
                if (test == NULL)
                  { fprintf(stderr,"%s: Cannot open %s/%s.1aln for output\n",
                                   Prog_Name,ONE_PATH,ONE_ROOT);
                    exit (1);
                  }
                fclose(test);
                break;
              }
            fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
            exit (1);
          case 'c':
            ARG_NON_NEGATIVE(CHAIN_MIN,"minimum seed cover");
            CHAIN_MIN <<= 1;
            break;
          case 'f':
            ARG_NON_NEGATIVE(FREQ,"maximum seed frequency");
            break;
          case 'i':
            ARG_REAL(ALIGN_RATE);
            if (ALIGN_RATE < .55 || ALIGN_RATE >= 1.)
              { fprintf(stderr,"%s: '-i' minimum alignment similarity must be in [0.55,1.0)\n",
                               Prog_Name);
                exit (1);
              }
            ALIGN_RATE = 1.-ALIGN_RATE;
            break;
          case 'l':
            ARG_NON_NEGATIVE(ALIGN_MIN,"minimum alignment length");
            break;
          case 'p':
            if (strncmp(argv[i]+1,"paf",3) == 0)
              { OUT_TYPE = 0;
                eptr = argv[i]+4;
                while (*eptr)
                  { switch (*eptr)
                      { case 'm':
                          OUT_OPT |= PAFM;
                          break;
                        case 'x':
                          OUT_OPT |= PAFX;
                          break;
                        case 's':
                          OUT_OPT |= PAFS;
                          break;
                        case 'S':
                          OUT_OPT |= PAFL;
                          break;
                        default:
                          fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
                          exit (1);
                          break;
                      }
                    eptr++;
                  }
                  break;
              }
            else if (strcmp(argv[i]+1,"psl") == 0)
              { OUT_TYPE = 1;
                break;
              }
            fprintf(stderr,"%s: Do not recognize option %s\n",Prog_Name,argv[i]);
            exit (1);
          case 's':
            ARG_NON_NEGATIVE(CHAIN_BREAK,"seed chain break threshold");
            CHAIN_BREAK <<= 1;
            break;
          case 'L':
            if (argv[i][2] != ':')
              { fprintf (stderr,"%s: option -L must be followed by :<filename>\n",Prog_Name);
                exit (1);
              }
            LOG_PATH = argv[i]+3;
            LOG_FILE = fopen (LOG_PATH, "a");
            if (LOG_FILE == NULL)
              { fprintf (stderr,"%s: Cannot open logfile %s for output\n",Prog_Name,LOG_PATH);
                exit(1);
              }
            break;
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
        }
      else
        argv[j++] = argv[i];
    argc = j;

    VERBOSE   = flags['v'];
    KEEP      = flags['k'];
    SOFT_MASK = flags['M'];
    SYMMETRIC = flags['S'];
```

---

## 3. 各选项详细说明

### 3.1 `-v` (Verbose) 选项

**处理位置：**
- **标志设置：** 第 4470 行 `ARG_FLAGS("vkMS")` - 将 `-v` 标记到 flags 数组
- **变量赋值：** 第 4567 行 `VERBOSE = flags['v'];`
- **变量声明：** 第 70 行 `static int VERBOSE;`

**功能：** 启用详细输出模式，显示处理过程的统计信息

---

### 3.2 `-k` (Keep) 选项

**处理位置：**
- **标志设置：** 第 4470 行 `ARG_FLAGS("vkMS")`
- **变量赋值：** 第 4568 行 `KEEP = flags['k'];`
- **变量声明：** 第 79 行 `static int KEEP;`

**功能：** 保留生成的 `.1gdb` 和 `.gix` 中间文件（默认会在完成后删除）

**使用位置：** 第 151-195 行的 `Clean_Exit()` 函数中检查 `KEEP` 标志决定是否删除中间文件

---

### 3.3 `-L:<log:path>` 选项

**处理位置：** 第 4544-4555 行

```4544:4555:FastGA.c
          case 'L':
            if (argv[i][2] != ':')
              { fprintf (stderr,"%s: option -L must be followed by :<filename>\n",Prog_Name);
                exit (1);
              }
            LOG_PATH = argv[i]+3;
            LOG_FILE = fopen (LOG_PATH, "a");
            if (LOG_FILE == NULL)
              { fprintf (stderr,"%s: Cannot open logfile %s for output\n",Prog_Name,LOG_PATH);
                exit(1);
              }
            break;
```

**变量声明：**
- 第 71 行 `static char *LOG_PATH;`
- 第 72 行 `static FILE *LOG_FILE;`

**功能：** 
- 指定日志文件路径（必须用 `:` 分隔，如 `-L:log.txt`）
- 以追加模式（`"a"`）打开日志文件
- 如果无法打开文件则退出程序

**使用位置：** 在整个程序中，当需要写入日志时检查 `LOG_FILE` 是否非空

---

### 3.4 `-P<dir>` 选项

**处理位置：** 第 4556-4558 行

```4556:4558:FastGA.c
          case 'P':
            SORT_PATH = argv[i]+2;
            break;
```

**变量声明：** 第 78 行 `static char *SORT_PATH;`

**默认值：** 第 4453-4455 行
```4453:4455:FastGA.c
    SORT_PATH = getenv("TMPDIR");
    if (SORT_PATH == NULL)
      SORT_PATH = ".";
```

**功能：** 
- 指定临时文件目录
- 如果未指定，使用环境变量 `TMPDIR`，如果 `TMPDIR` 未设置则使用当前目录 `.`

**使用位置：** 
- 传递给 `GIXmake` 和 `FAtoGDB` 子程序（第 4748-4752, 4786-4790 行）
- 用于存储排序和索引过程中的临时文件

---

### 3.5 `-1:<align:path>` 选项

**处理位置：** 第 4472-4484 行

```4472:4484:FastGA.c
          case '1':
            if (strncmp(argv[i]+1,"1:",2) == 0)
              { OUT_TYPE = 2;
                ONE_PATH = PathTo(argv[i]+3);
                ONE_ROOT = Root(argv[i]+3,".1aln");
                test = fopen(Catenate(ONE_PATH,"/",ONE_ROOT,".1aln"),"w");
                if (test == NULL)
                  { fprintf(stderr,"%s: Cannot open %s/%s.1aln for output\n",
                                   Prog_Name,ONE_PATH,ONE_ROOT);
                    exit (1);
                  }
                fclose(test);
                break;
              }
```

**变量声明：**
- 第 85 行 `static char *ONE_PATH;`
- 第 86 行 `static char *ONE_ROOT;`
- 第 83 行 `static int OUT_TYPE;` （设置为 2 表示 ONEcode ALN 格式）

**功能：**
- 指定输出 ONEcode ALN 格式文件
- 格式必须为 `-1:路径`（冒号分隔）
- 自动添加 `.1aln` 扩展名（如果未提供）
- 在解析时测试文件是否可以创建

**使用位置：** 在写入对齐结果时使用（第 4922-4927 行）

---

### 3.6 `-T<int>` 选项

**处理位置：** 第 4559-4561 行

```4559:4561:FastGA.c
          case 'T':
            ARG_NON_NEGATIVE(NTHREADS,"number of threads to use");
            break;
```

**变量声明：** 第 77 行 `static int NTHREADS;`

**默认值：** 第 4456 行 `NTHREADS = 8;`

**功能：** 指定使用的线程数

---

## 4. 帮助信息输出

**位置：** 第 4572-4599 行

当参数数量不正确时，程序会输出使用说明：

```4572:4599:FastGA.c
    if (argc != 3 && argc != 2)
      { fprintf(stderr,"\nUsage: %s %s\n",Prog_Name,Usage[0]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[1]);
        fprintf(stderr,"       %*s %s\n",(int) strlen(Prog_Name),"",Usage[2]);
        fprintf(stderr,"\n");
        fprintf(stderr,"         <format> = -paf[mxsS]* | -psl | -1:<align:path>[.1aln]\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"         <precursor> = .gix | .1gdb | <fa_extn> | <1_extn>\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"             <fa_extn> = (.fa|.fna|.fasta)[.gz]\n");
        fprintf(stderr,"             <1_extn>  = any valid 1-code sequence file type\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -v: Verbose mode, output statistics as proceed.\n");
        fprintf(stderr,"      -k: Keep any generated .1gdb's and .gix's.\n");
        fprintf(stderr,"      -M: Use soft mask information if available.\n");
        fprintf(stderr,"      -S: Use symmetric seeding (not recommended).\n");
        fprintf(stderr,"      -L: Output log to specified file.\n");
        fprintf(stderr,"      -T: Number of threads to use.\n");
        fprintf(stderr,"      -P: Directory to use for temporary files.\n");
        fprintf(stderr,"\n");
        fprintf(stderr,"      -paf: Stream PAF output\n");
        fprintf(stderr,"        -pafx: Stream PAF output with CIGAR string with X's\n");
        fprintf(stderr,"        -pafm: Stream PAF output with CIGAR string with ='s\n");
        fprintf(stderr,"        -pafs: Stream PAF output with CS string in short form\n");
        fprintf(stderr,"        -pafS: Stream PAF output with CS string in long form\n");
        fprintf(stderr,"      -psl: Stream PSL output\n");
        fprintf(stderr,"      -1: Generate 1-code output to specified file\n");
        fprintf(stderr,"\n");
```

---

## 5. 辅助函数

### ARG_FLAGS 宏
用于处理简单的标志选项（如 `-v`, `-k`, `-M`, `-S`），这些选项不需要参数值。

### ARG_NON_NEGATIVE 宏
用于解析非负整数参数（如 `-T8`, `-f10`）。

### PathTo() 和 Root() 函数
用于从完整路径中提取目录路径和文件名根。

---

## 总结

| 选项 | 代码位置 | 关键行号 | 功能 |
|------|---------|---------|------|
| `-v` | FastGA.c | 4470, 4567 | 详细输出模式 |
| `-k` | FastGA.c | 4470, 4568 | 保留中间文件 |
| `-L:<path>` | FastGA.c | 4544-4555 | 指定日志文件 |
| `-P<dir>` | FastGA.c | 4556-4558 | 指定临时文件目录 |
| `-1:<path>` | FastGA.c | 4472-4484 | 指定ALN输出文件 |
| `-T<int>` | FastGA.c | 4559-4561 | 指定线程数 |

所有选项的解析都在 `main()` 函数的参数处理循环中（第 4466-4565 行）完成。

