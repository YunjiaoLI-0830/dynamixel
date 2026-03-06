import threading
import queue
import tkinter as tk
from tkinter import ttk

from dynamixel_controller import XM430Controller, DynamixelConfig, DynamxielError

class App(tk.Tk):
    """
    一个简单的 Tkinter GUI，用于控制 XM430 线性机构：
    - Connect / Disconnect（连接串口）
    - Set Home（把当前位置记为统一起点 HOME）
    - Reset to Home（回到 HOME）
    - Move（输入“相对 HOME 的位移(mm)”并移动）
    - Read Position（读取当前 ticks）
    重点：UI 与电机控制分离，电机控制在 dynamixel_controller.py 中。
    """

    def __init__(self):
        super().__init__()
        self.title("XM430 Linear Actuator Controller")
        self.geometry("760x520")

        # 消息队列：后台线程把日志/错误放进来，主线程定时取出显示到文本框
        # 目的：避免线程直接操作 Tkinter（Tkinter 只能在主线程更新 UI）
        self.msg_q: queue.Queue[str] = queue.Queue()

        # 控制器实例（控制层）：负责串口通信、home、move 等
        # 默认参数在 DynamixelConfig() 中，UI 会允许你修改
        self.ctrl = XM430Controller(DynamixelConfig())

        # ---------- UI 变量（tkinter 的变量类） ----------
        # 优点：Entry 与变量绑定，读取/写入会自动同步
        self.port_var = tk.StringVar(value=self.ctrl.cfg.port)
        self.baud_var = tk.IntVar(value=self.ctrl.cfg.baudrate)
        self.id_var = tk.IntVar(value=self.ctrl.cfg.dxl_id)

        self.lead_var = tk.DoubleVar(value=self.ctrl.cfg.lead_mm_per_rev)
        self.gear_var = tk.DoubleVar(value=self.ctrl.cfg.gear_ratio)
        self.maxtravel_var = tk.DoubleVar(value=self.ctrl.cfg.max_travel_mm)
        self.tol_var = tk.DoubleVar(value=self.ctrl.cfg.inpos_tol_mm)
        self.timeout_var = tk.DoubleVar(value=self.ctrl.cfg.timeout_s)

        # 运动输入：相对 HOME 的位移（mm）
        self.move_mm_var = tk.DoubleVar(value=0.4)

        # 搭建界面组件
        self._build_ui()

        # 开始轮询消息队列，把后台线程输出写进 log 文本框
        self._poll_queue()

        # 关闭窗口时：安全断开电机 & 退出程序
        self.protocol("WM_DELETE_WINDOW", self.on_close)

    def _build_ui(self):
        """搭建 UI 布局：Connection / Params / Motion / Log 四块"""

        frm = ttk.Frame(self, padding=10)
        frm.pack(fill="both", expand=True)

        # ----------------- 1) Connection 区（串口连接） -----------------
        g1 = ttk.LabelFrame(frm, text="Connection", padding=10)
        g1.pack(fill="x")

        ttk.Label(g1, text="Port").grid(row=0, column=0, sticky="w")
        ttk.Entry(g1, textvariable=self.port_var, width=12).grid(row=0, column=1, padx=5)

        ttk.Label(g1, text="Baud").grid(row=0, column=2, sticky="w")
        ttk.Entry(g1, textvariable=self.baud_var, width=10).grid(row=0, column=3, padx=5)

        ttk.Label(g1, text="ID").grid(row=0, column=4, sticky="w")
        ttk.Entry(g1, textvariable=self.id_var, width=6).grid(row=0, column=5, padx=5)

        ttk.Button(g1, text="Connect", command=self.connect).grid(row=0, column=6, padx=10)
        ttk.Button(g1, text="Disconnect", command=self.disconnect).grid(row=0, column=7)

        # ----------------- 2) Params 区（机构参数/软限位/阈值） -----------------
        g2 = ttk.LabelFrame(frm, text="Mechanism Params", padding=10)
        g2.pack(fill="x", pady=10)

        row = 0
        ttk.Label(g2, text="Lead (mm/rev)").grid(row=row, column=0, sticky="w")
        ttk.Entry(g2, textvariable=self.lead_var, width=10).grid(row=row, column=1, padx=5)

        ttk.Label(g2, text="Gear Ratio").grid(row=row, column=2, sticky="w")
        ttk.Entry(g2, textvariable=self.gear_var, width=10).grid(row=row, column=3, padx=5)

        ttk.Label(g2, text="Max Travel (mm)").grid(row=row, column=4, sticky="w")
        ttk.Entry(g2, textvariable=self.maxtravel_var, width=10).grid(row=row, column=5, padx=5)

        row = 1
        ttk.Label(g2, text="InPos Tol (mm)").grid(row=row, column=0, sticky="w")
        ttk.Entry(g2, textvariable=self.tol_var, width=10).grid(row=row, column=1, padx=5)

        ttk.Label(g2, text="Timeout (s)").grid(row=row, column=2, sticky="w")
        ttk.Entry(g2, textvariable=self.timeout_var, width=10).grid(row=row, column=3, padx=5)

        # Apply Params：把 UI 上输入的参数写回 controller.cfg（控制层）
        ttk.Button(g2, text="Apply Params", command=self.apply_params).grid(row=1, column=5, padx=5)

        # ----------------- 3) Motion 区（HOME-based 运动） -----------------
        g3 = ttk.LabelFrame(frm, text="Motion (HOME-based)", padding=10)
        g3.pack(fill="x")

        # Set Home：把当前 ticks 记为 HOME（统一起始点）
        ttk.Button(g3, text="Set Home", command=self.set_home).grid(row=0, column=0, padx=5)

        # Reset to Home：回到 HOME
        ttk.Button(g3, text="Reset to Home", command=self.reset_home).grid(row=0, column=1, padx=5)

        ttk.Label(g3, text="Move to (mm from Home)").grid(row=0, column=2, padx=10, sticky="e")

        # 输入框：位移（相对 HOME 的绝对位置）
        ent = ttk.Entry(g3, textvariable=self.move_mm_var, width=10)
        ent.grid(row=0, column=3, padx=5)

        # 绑定回车键：按 Enter 就触发 Move
        ent.bind("<Return>", lambda e: self.move_from_home())

        ttk.Button(g3, text="Move", command=self.move_from_home).grid(row=0, column=4, padx=5)

        # Read Position：读当前 ticks（用于调试）
        ttk.Button(g3, text="Get Current Position", command=self.get_current_position).grid(row=0, column=5, padx=5)

        # ----------------- 4) Log / Errors 区（显示日志和错误） -----------------
        g4 = ttk.LabelFrame(frm, text="Log / Errors", padding=10)
        g4.pack(fill="both", expand=True, pady=10)

        self.log = tk.Text(g4, height=12, wrap="word")
        self.log.pack(fill="both", expand=True)
        self._log_line("Ready.")

        # 小提示
        hint = ttk.Label(frm, text="Tip: Connect 后先 Set Home；输入框按 Enter 可直接 Move。", foreground="gray")
        hint.pack(anchor="w")

    # ----------------- 日志相关工具函数 -----------------
    def _log_line(self, s: str):
        """在文本框追加一行日志"""
        self.log.insert("end", s + "\n")
        self.log.see("end")  # 自动滚动到底部

    def _run_bg(self, fn, *args):
        """
        关键：后台线程执行控制层操作，避免 UI 卡死。
        因为电机操作（open/move/wait）可能需要几十 ms~几秒。
        """
        def worker():
            try:
                fn(*args)
            except DynamxielError as e:
                # 控制层抛出的“可读”错误：直接显示在 UI
                self.msg_q.put(f"[ERROR] {e}")
            except Exception as e:
                # 其它异常（程序 bug 等）
                self.msg_q.put(f"[ERROR] {type(e).__name__}: {e}")

        threading.Thread(target=worker, daemon=True).start()

    def _poll_queue(self):
        """
        UI 主线程定时（每 100ms）检查队列是否有消息，有就写到 log。
        这样后台线程不会直接碰 UI，也不会崩。
        """
        while True:
            try:
                msg = self.msg_q.get_nowait()
            except queue.Empty:
                break
            self._log_line(msg)

        self.after(100, self._poll_queue)

    # ----------------- 把 UI 参数写回 controller.cfg -----------------
    def apply_params(self):
        """
        把 UI 输入的参数写回控制器配置。
        注意：如果已经连接了，改 port/baud 需要重新 connect 才生效。
        """
        self.ctrl.cfg.port = self.port_var.get().strip()
        self.ctrl.cfg.baudrate = int(self.baud_var.get())
        self.ctrl.cfg.dxl_id = int(self.id_var.get())

        self.ctrl.cfg.lead_mm_per_rev = float(self.lead_var.get())
        self.ctrl.cfg.gear_ratio = float(self.gear_var.get())
        self.ctrl.cfg.max_travel_mm = float(self.maxtravel_var.get())
        self.ctrl.cfg.inpos_tol_mm = float(self.tol_var.get())
        self.ctrl.cfg.timeout_s = float(self.timeout_var.get())

        self._log_line("[INFO] Params applied.")

    # ----------------- UI 按钮：Connect/Disconnect -----------------
    def connect(self):
        """Connect 按钮：更新参数 -> 后台 open()"""
        self.apply_params()
        self._log_line("[INFO] Connecting...")
        self._run_bg(self._connect_impl)

    def _connect_impl(self):
        self.ctrl.open()
        self.msg_q.put("[INFO] Connected. (建议先 Set Home)")

    def disconnect(self):
        """Disconnect 按钮：后台 close()"""
        self._log_line("[INFO] Disconnecting...")
        self._run_bg(self._disconnect_impl)

    def _disconnect_impl(self):
        self.ctrl.close()
        self.msg_q.put("[INFO] Disconnected.")

    # ----------------- UI 按钮：Home/Reset/Move -----------------
    def set_home(self):
        """Set Home 按钮：后台 set_home()"""
        self._run_bg(self._set_home_impl)

    def _set_home_impl(self):
        ticks = self.ctrl.set_home()
        self.msg_q.put(f"[HOME] set to {ticks} ticks")

    def reset_home(self):
        """Reset to Home 按钮：后台 reset_to_home()"""
        self._run_bg(self._reset_home_impl)

    def _reset_home_impl(self):
        self.ctrl.reset_to_home(wait=True)
        self.msg_q.put("[HOME] reset done")

    def move_from_home(self):
        """Move 按钮 / Enter：读取输入 mm -> 后台 move_from_home_mm()"""
        mm = float(self.move_mm_var.get())
        self._run_bg(self._move_from_home_impl, mm)

    def _move_from_home_impl(self, mm):
        self.msg_q.put(f"[MOVE] target = HOME {mm:+.3f} mm")
        self.ctrl.move_from_home_mm(mm, wait=True)
        self.msg_q.put("[MOVE] done")

    # >>> ADDED / MODIFIED: 新函数名与更详细输出
    def get_current_position(self):
        self._run_bg(self._get_current_position_impl)

    def _get_current_position_impl(self):
        ticks = self.ctrl.get_present_ticks()

        # 如果已经 set home，则显示“相对home的mm”
        if self.ctrl.home_ticks is not None:
            rel_ticks = ticks - self.ctrl.home_ticks
            # mm 反算：mm = ticks / (ticks_per_rev*gear_ratio) * lead
            mm = (rel_ticks / (self.ctrl.cfg.ticks_per_rev * self.ctrl.cfg.gear_ratio)) * self.ctrl.cfg.lead_mm_per_rev
            self.msg_q.put(f"[POS] present = {ticks} ticks | from HOME = {mm:+.3f} mm")
        else:
            self.msg_q.put(f"[POS] present = {ticks} ticks (HOME not set)")

    # ----------------- 关闭窗口时的清理 -----------------
    def on_close(self):
        """
        关闭窗口：尝试断开电机（Torque Off + close port），再退出 UI。
        """
        try:
            self.ctrl.close()
        except Exception:
            pass
        self.destroy()


if __name__ == "__main__":
    App().mainloop()