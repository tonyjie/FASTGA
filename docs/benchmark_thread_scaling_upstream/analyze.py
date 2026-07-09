import csv, statistics as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN="/tmp/claude-1564080/-work-shared-users-phd-jl4257-Project-genomics-agent-FASTGA/66f02e18-e698-4955-9f60-687e7a3925eb/scratchpad/tscaling"

rows=[]
with open(f"{RUN}/results.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["nonredundant_aln"] and int(r["threads"])<=32:
            rows.append(r)

by_t={}
for r in rows:
    t=int(r["threads"]); by_t.setdefault(t,[]).append(r)

T=sorted(by_t)
def med(t,k): return st.median(float(r[k]) for r in by_t[t])
wall={t:med(t,"wall_s") for t in T}
cpu ={t:med(t,"cpu_pct") for t in T}
rss ={t:med(t,"maxrss_kb")/1024 for t in T}   # MB
base=wall[1]
speed={t:base/wall[t] for t in T}

# old build (optimize-memory) medians, seconds
old={1:133.68,2:73.13,4:48.83,8:32.63,16:24.60,32:22.03}
old_speed={t:old[1]/old[t] for t in old}

# ---- print markdown table ----
print("| T | wall(s) new | speedup new | CPU% | RSS(MB) | wall(s) old | speedup old |")
print("|--:|--:|--:|--:|--:|--:|--:|")
for t in T:
    print(f"| {t} | {wall[t]:.1f} | {speed[t]:.2f}x | {cpu[t]:.0f} | {rss[t]:.0f} | {old[t]:.1f} | {old_speed[t]:.2f}x |")

# ---- plot ----
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.8))
C_NEW="#2563eb"; C_OLD="#9ca3af"; C_IDEAL="#d1d5db"

ax1.plot(T,[wall[t] for t in T],"o-",color=C_NEW,lw=2,label="upstream 10ebff7 (new)")
ax1.plot(list(old),[old[t] for t in old],"s--",color=C_OLD,lw=1.6,label="optimize-memory (old)")
ax1.set_xscale("log",base=2); ax1.set_yscale("log",base=2)
ax1.set_xticks(T); ax1.set_xticklabels(T)
ax1.set_xlabel("threads (-T)"); ax1.set_ylabel("wall-clock (s, log2)")
ax1.set_title("FastGA end-to-end runtime\nEXAMPLE HAP1 vs HAP2 (~86 Mbp each)")
ax1.grid(True,which="both",alpha=0.25); ax1.legend(frameon=False)

ax2.plot(T,[speed[t] for t in T],"o-",color=C_NEW,lw=2,label="new (measured)")
ax2.plot(list(old),[old_speed[t] for t in old],"s--",color=C_OLD,lw=1.6,label="old (measured)")
ax2.plot(T,T,":",color=C_IDEAL,lw=1.6,label="ideal linear")
ax2.set_xscale("log",base=2); ax2.set_yscale("log",base=2)
ax2.set_xticks(T); ax2.set_xticklabels(T); ax2.set_yticks(T); ax2.set_yticklabels(T)
ax2.set_xlabel("threads (-T)"); ax2.set_ylabel("speedup vs T1 (log2)")
ax2.set_title("Parallel speedup (32-thread hard cap)")
ax2.grid(True,which="both",alpha=0.25); ax2.legend(frameon=False)

fig.tight_layout()
out=f"{RUN}/thread_scaling_new_upstream.png"
fig.savefig(out,dpi=140,bbox_inches="tight")
print("\nsaved:",out)
