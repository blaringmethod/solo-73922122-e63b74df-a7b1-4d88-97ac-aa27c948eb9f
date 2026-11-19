# VCF 队列合并服务（FastAPI + SQLite + pysam）

把分批检测得到的 2~3 份二倍体 VCF 合并成一份队列 VCF，供关联分析使用。
网页完成上传、校验、重复样本取用源选择、后台生成与结果下载，全程在
SQLite 中记录状态，并通过“选择版本号”保证下载到的永远是**当前选择对应
的完整结果**。

## 一、业务规则实现对照

| 需求 | 实现方式 |
| --- | --- |
| 仅接受 2~3 份 VCF（.vcf / .vcf.gz） | `POST /api/batches` 数量与扩展名校验 |
| 同参考版本、坐标排序、每份同坐标仅一条记录 | `app/validation.py` 单文件校验 + `app/cross_check.py` 跨文件检查 |
| 等位基因限 ACGT 序列 | REF/ALT 正则 `^[ACGT]+$`（ALT 为 `.` 的无 ALT 位点合法） |
| FORMAT 限 GT、AD、PL，且 GT 必需 | 逐记录检查 FORMAT 字段集合 |
| 二倍体 | 每个样本 GT 必须是两个等位槽位（允许 `./.` 缺失） |
| 重复样本必须选择唯一取用原件，未选定不能生成 | `dup_groups` 表 + 生成前置校验 `assert_generation_allowed` |
| 整份检测不能逐位点拼接 | 每个样本在整份输出中只来自一个被选原件（`SampleSource`），缺失位点保持缺失 |
| 联合各位点 ALT，同步重排 GT/AD/PL | `app/merge.py`：ALT 按首次出现取并集，GT 索引、AD（按等位）、PL（三角布局）整体重映射 |
| 新增等位的未知值保留缺失 | 其他文件引入的等位对应 AD/PL 槽位输出 `.`，不伪造数值 |
| 缺少位点不能当作 0/0 | 源文件无该坐标时输出 `./.` 且 AD/PL 全缺失 |
| 同坐标 REF 冲突则整份失败 | 跨文件 REF 比对，冲突写入批次级错误并阻断生成 |
| 换件只能替换同批一份原件，其余原件须已通过单文件校验 | `POST /api/batches/{id}/replacements` 同步检查目标属于该批、其余原件 `valid=1` |
| 换件须校验新原件及其换件后的整批文件 | 后台校验新原件单文件规则，再以“新原件 + 其余原件”跑跨文件检查（参考版本、REF 冲突、contig 顺序） |
| 校验期间仍可改选样本来源 | 选择接口不因换件校验而禁用；提交换件的事务内重读最新选择 |
| 校验期间换件阻止生成 | 校验中 `ready_to_generate=false`，生成接口返回 409 |
| 校验失败显示原因，原批次照常使用 | 换件记录 `failed` 并带原因；原件、校验状态、分组、版本号均不变 |
| 校验期间有另一次换件先提交，本次拒绝并提示刷新 | 每次换件记录提交时的 `replace_seq`；提交事务内不一致即 `rejected`，页面提示刷新 |
| 成功换件一并更新原件、校验状态与样本分组 | `commit_replacement` 单事务内落库并重建 `dup_groups` |
| 旧合并结果及在途生成立即失效 | 同一事务内 `selection_version+1`，done/running/pending 全部 `superseded`，输出文件删除 |
| 换件后只有一个来源的样本自动选定 | 重建分组时 `len(file_ids)==1` 直接选定 |
| 换件前多来源且已选中未换原件则沿用，否则待确认 | 重建分组时仅当“换件前已多来源”且“原选择仍是未换原件的来源”才沿用，否则清空选择 |
| 后台生成时可改取用源，改动立即使旧结果及在途生成失效 | 每次变更 `selection_version+1`，既有 generation 立即置 `superseded`，worker 在位点间检查版本并中止 |
| 改回也不能恢复 | 版本号只增不减；旧 generation 不会回退，输出文件被删除，需重新生成 |
| 只下载当前选择对应的完整 VCF | 下载接口校验最新 generation 为 done、版本号等于当前版本、文件存在 |
| 失败或重启不开放半成品 | 合并先写 `.partial`，成功才原子 `os.replace`；异常/中止即删除；启动时把在途任务标记失败 |
| 网页显示样本来源和失败原因 | 单页应用展示每份文件的样本、位点、状态、失败原因及批次级错误 |

## 二、目录结构

```
.
├── app/
│   ├── config.py        # 路径、常量、状态枚举
│   ├── db.py            # SQLite schema、连接、版本失效、换件提交事务、启动恢复
│   ├── validation.py    # 单文件 VCF 校验
│   ├── cross_check.py   # 参考版本 / REF 冲突 / contig 顺序
│   ├── merge.py         # ALT 联合与 GT/AD/PL 重排合并引擎
│   ├── service.py       # 上传编排、选择、换件后台校验、后台 worker
│   ├── main.py          # FastAPI 路由与 JSON 序列化
│   └── static/index.html# 单页前端
├── demo/
│   ├── generate_demo.py # 演示文件生成脚本
│   ├── batch_a.vcf      # 样本 S1,S2（多等位、无 ALT、缺失基因型）
│   ├── batch_b.vcf      # 样本 S2,S3（ALT 顺序不同、独有位点）
│   ├── batch_c.vcf      # 样本 S3,S4
│   ├── replacement_b.vcf    # batch_b 的更正补发件（S2,S5），用于换件演示
│   ├── bad_replacement_b.vcf # 与 batch_a 在 chr1:100 REF 冲突，用于换件失败演示
│   └── bad_ref_conflict.vcf # 故意制造 chr1:100 REF 冲突的反例
├── tests/               # 42 个 pytest 用例
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── run.py               # 本地启动入口
```

## 三、本地运行

需要 Python 3.11（依赖见 `requirements.txt`）：

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python run.py
# 打开 http://localhost:8000
```

运行测试（先安装开发依赖 `pip install -r requirements-dev.txt`）：

```bash
.venv/bin/python -m pytest tests/ -q
```

重新生成演示文件：

```bash
.venv/bin/python demo/generate_demo.py demo
```

## 四、Docker Compose 启动（可实际启动）

```bash
docker compose up --build
# 浏览器打开 http://localhost:8000
```

- 服务监听容器内 8000，映射到宿主机 8000；
- SQLite 数据库、上传原件、输出 VCF 与临时文件都在命名卷 `vcf_data`
  （容器内 `/data`）下，重启容器不丢批次；
- 容器健康检查调用 `/api/health`。

> 说明：本环境未安装 Docker 守护进程，Dockerfile/Compose 未在本机实际
> 构建，但镜像基于官方 `python:3.11-slim`，pysam 使用 manylinux 预编译
> wheel，无需系统级编译依赖。

## 五、操作流程

1. 网页选择 2~3 份 VCF，点“上传并校验”。
2. 查看每份文件的状态、样本、位点数与失败原因；批次级错误（参考版本
   不一致、REF 冲突等）在页面顶部红色区域展示。
3. 在“重复样本取用源选择”中为每个出现在多份文件中的样本指定唯一原件；
   只有一个来源的样本自动选定。所有重复样本选定前“开始生成”按钮禁用。
4. 点“开始生成”，页面每 1.5 秒轮询状态（等待中/生成中/已完成/失败/已失效）。
5. 完成后点“下载当前选择对应的完整 VCF”。
6. 任何时候改动取用源：页面提示旧结果与在途任务已失效，旧输出文件立即
   删除；需要重新生成后才能再次下载。

### 换件演示（实验室补发更正原件）

1. 上传 `demo/batch_a.vcf`、`demo/batch_b.vcf`（S2 同时出现在两份文件中）。
2. 在“4. 换件”中选择要替换的原件（例如 `batch_b.vcf`），上传补发的
   `demo/replacement_b.vcf`，点“提交换件并校验”。
3. 校验期间（Compose 默认人为延迟 4 秒，见
   `VCF_REPLACE_VALIDATE_DELAY_SECONDS`）仍可改“取用源选择”。
4. 结果：
   - 成功时该原件被替换，`S2` 若原先选中的是未换原件则沿用，否则变为
     “待确认”；只有单一来源的样本自动选定；旧结果标记“已失效”。
   - 上传 `demo/bad_replacement_b.vcf`（REF 与 A 冲突）则换件失败并显示
     原因，原批次仍可照常使用。
   - 校验期间连续提交两次换件，先完成者提交，后者显示“已被拒绝，请刷新”。

## 六、HTTP API 摘要

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/batches` | multipart 上传 2~3 份 VCF，返回批次全量状态 |
| GET | `/api/batches/{id}` | 查询批次、文件、选择分组、换件记录、最新 generation |
| POST | `/api/batches/{id}/replacements` | multipart 字段 `target_file_id` + `file`，后台校验换件并返回批次状态（含 `replacement_id`） |
| PUT | `/api/batches/{id}/selection/{group_id}` | 表单字段 `file_id`（空串表示取消选择），变更即 bump 版本并失效旧任务 |
| POST | `/api/batches/{id}/generate` | 前置校验通过后启动后台合并线程 |
| GET | `/api/batches/{id}/download` | 仅当最新任务 done 且版本匹配时返回完整 VCF |
| GET | `/api/health` | 健康检查 |

## 七、合并语义细节

- 输出 header 基于第一份源文件重建（contig、INFO/FORMAT/FILTER 定义），
  样本列替换为去重后的队列样本（按源文件顺序首次出现排序），并追加
  `##source` 记录。
- 同一坐标的 ALT 取并集，顺序按文件顺序、位点内 ALT 顺序的首次出现。
- GT：源等位索引映射到联合索引；保留原相位（`/` 与 `|`）。
- AD：按等位逐槽位移位，联合中由其他文件引入的等位槽位置 `.`。
- PL：按 VCF 二倍体三角顺序 `(00),(01),(11),(02),(12),(22)...`
  （索引 `j*(j+1)/2+i`）重映射，未知基因型位置 `.`。
- 源文件缺少该坐标的样本：`GT=./.`、AD/PL 全缺失，绝不写 `0/0`。
- 输出固定 `FORMAT=GT:AD:PL`，QUAL/FILTER/INFO 留空（仅保留 header 定义）。

## 八、测试覆盖（42 项）

- 单文件校验 10 项：多等位合法文件、未排序、重复坐标、非 ACGT、非法
  FORMAT、缺 GT、非二倍体、缺失二倍体基因型、空文件、未声明 contig；
- 合并引擎 8 项：ALT 联合与 GT/AD/PL 重排、新等位缺失值、缺失位点不
  变 0/0、重复样本只取一个原件、相位保留、无 ALT 位点、中止清理、PL
  索引布局；
- API/工作流 11 项：健康检查、上传数量、完整演示流程、REF 冲突整批失败、
  单文件损坏原因、选择变更失效已完成结果（含改回不恢复）、在途任务被
  中止且无半成品、重启恢复阻断下载、样本来源展示、BGZF 压缩 VCF 上传、
  非法扩展名拒绝；
- 换件 13 项（`tests/test_replacement.py`）：
  - 端到端 7 项：成功换件同时更新原件/校验状态/分组、被换原件恰为已选
    来源时改为待确认、换件失效已完成结果并删除旧输出、换件后整批校验
    失败时显示原因且原批次照常可用、新原件单文件校验失败、其余原件未
    通过校验时同步拒绝、校验期间可改选样本来源且提交时按最新选择重建；
  - 并发/恢复 2 项：校验期间另一次换件先提交则本次被拒绝（刷新提示）、
    重启把在途换件判定失败；
  - 分组重建规则 4 项：多来源且选中未换原件则沿用、选中被换原件则清空、
    换件前单来源换件后多来源须待确认、单来源自动选定。

## 九、关于 .vcf.gz

pysam/htslib 要求 `.vcf.gz` 为 **BGZF**（`bgzip`）格式而非普通 gzip，
否则无法随机访问并报 `seek not implemented`。请使用：

```bash
bgzip -c input.vcf > input.vcf.gz
```

仓库 `tests/bgzf_util.py` 是一个仅依赖标准库 zlib 的 BGZF 写入器，供测试
与演示生成合规 `.vcf.gz`，生产环境仍建议使用 htslib 自带的 `bgzip`。
上传普通 gzip 时，页面会在该文件的失败原因中明确提示需要 BGZF。
