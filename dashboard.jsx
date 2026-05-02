import React, { useState, useEffect, useMemo, useCallback } from "react";
import { LineChart, Line, BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, PieChart, Pie, Cell } from "recharts";
import _ from "lodash";

// ── Constantes para los menús desplegables ──
const SUCURSALES = ["Bernardo Quintana", "Belen", "Centro sur"];
const EMPLEADOS = ["Ana García", "Carlos López", "María Sánchez", "Pedro Ramírez", "Laura Torres", "Diego Hernández"];
const PRODUCTOS = ["LENTES SOL", "LENTES GRADUADOS", "ARMAZÓN", "MICA", "ESTUCHE", "SOLUCIÓN", "LENTES CONTACTO", "ACCESORIOS", "BIFOCAL", "PROGRESIVO", "ANTIREFLEJANTE", "FILTRO AZUL", "CADENA", "MICROFIBRA", "GOTAS"];
const RAZONES_NO = ["PRECIO", "NO ENCONTRÓ LO QUE BUSCABA", "SOLO COTIZANDO", "REGRESA DESPUÉS", "NO LE GUSTÓ", "MUY CARO", "SIN GRADUACIÓN", "COMPARANDO PRECIOS", "FALTA DE STOCK", "NO TRAÍA RECETA"];
const COMENTARIOS_LIST = ["BUEN SERVICIO", "VARIEDAD LIMITADA", "PRECIOS ACCESIBLES", "REGRESA PRONTO", "EXCELENTE ATENCIÓN", "FALTA SURTIDO", "BUENA EXPERIENCIA", "MUY RÁPIDO", "RECOMENDARÍA", "NECESITA MÁS OPCIONES"];
const SEXOS = ["H", "M"];
const RANGOS_EDAD = ["Menor", "18-25", "26-35", "36-50", "51-65", "65+"];
const PAGOS = ["EFE", "TC"];
const MONTH_NAMES = ["Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"];

// ── AQUÍ EMPIEZA TU DASHBOARD REAL ──
export default function SalesDashboard() {
  // 1. Estado para guardar los datos de PostgreSQL
  const [allData, setAllData] = useState([]);
  const [loading, setLoading] = useState(true);

  // 2. Traer los datos de tu API (Node.js) al cargar la página
  useEffect(() => {
    const fetchRealData = async () => {
      try {
        const response = await fetch('http://localhost:3001/api/datos');
        const realData = await response.json();
        
        // Aseguramos el formato de fecha
        const formattedData = realData.map(item => ({
          ...item,
          fecha: new Date(item.fecha).toISOString().split("T")[0]
        }));

        setAllData(formattedData);
        setLoading(false);
      } catch (error) {
        console.error("Error cargando datos de PostgreSQL:", error);
        setLoading(false);
      }
    };

    fetchRealData();
  }, []);

  // 3. Pantalla de carga mientras trae los datos
  if (loading) {
    return (
      <div style={{ minHeight: "100vh", background: "#0c0c1d", color: "#e2e8f0", display: "flex", justifyContent: "center", alignItems: "center", fontSize: "20px" }}>
        Cargando datos desde Railway... ⏳
      </div>
    );
  }

  // 👇 DE AQUÍ PARA ABAJO, DEJA TUS ESTADOS COMO LOS TENÍAS 👇
  const [activeTab, setActiveTab] = useState("resumen");
  // const [sucursal, setSucursal] = useState("Todas");
  const [sucursal, setSucursal] = useState("Todas");
  const [kpi, setKpi] = useState("monto");
  const isMoney = kpi === "monto";

  // Date filters
  const [dateMode, setDateMode] = useState("todo");
  const [selectedYear, setSelectedYear] = useState(2026);
  const [selectedMonth, setSelectedMonth] = useState(new Date().getMonth());
  const [exactDate, setExactDate] = useState("");
  const [rangeStart, setRangeStart] = useState("");
  const [rangeEnd, setRangeEnd] = useState("");

  // Advanced filters
  const [showFilters, setShowFilters] = useState(false);
  const [filterSexo, setFilterSexo] = useState([]);
  const [filterEdad, setFilterEdad] = useState([]);
  const [filterPago, setFilterPago] = useState([]);
  const [filterFrecuente, setFilterFrecuente] = useState([]);
  const [filterEmpleado, setFilterEmpleado] = useState([]);
  const [filterPrecioMin, setFilterPrecioMin] = useState(0);
  const [filterPrecioMax, setFilterPrecioMax] = useState(5000);

  // Visible charts
  const [visibleCharts, setVisibleCharts] = useState(["ventas_hora", "ventas_mes", "ventas_semana", "top_productos"]);

  // Comparison layers
  const [horaLayers, setHoraLayers] = useState([]);
  const [mesLayers, setMesLayers] = useState([]);
  const [compMesYear, setCompMesYear] = useState(2026);
  const [compMesMonth, setCompMesMonth] = useState(new Date().getMonth());

  // Table tab
  const [tableSucursal, setTableSucursal] = useState("Todas");
  const [tableCompro, setTableCompro] = useState(["SI", "NO"]);
  const [tablePage, setTablePage] = useState(0);
  const pageSize = 30;

  // ── Filter data ──
  const filteredData = useMemo(() => {
    let d = [...allData]; // <--- AQUI CAMBIAMOS ALL_DATA por allData
    if (sucursal !== "Todas") d = d.filter(r => r.sucursal === sucursal);
    if (dateMode === "mes") d = d.filter(r => r.ano === Number(selectedYear) && r.mes_num === selectedMonth + 1);
    else if (dateMode === "fecha" && exactDate) d = d.filter(r => r.fecha === exactDate);
    else if (dateMode === "rango" && rangeStart && rangeEnd) d = d.filter(r => r.fecha >= rangeStart && r.fecha <= rangeEnd);
    if (filterSexo.length) d = d.filter(r => filterSexo.includes(r.sexo_h_m) || (filterSexo.includes("No especificado") && !r.sexo_h_m));
    if (filterEdad.length) d = d.filter(r => filterEdad.includes(r.rango_edad) || (filterEdad.includes("No especificado") && !r.rango_edad));
    if (filterPago.length) d = d.filter(r => filterPago.includes(r.pago_tc_efe) || (filterPago.includes("No especificado") && !r.pago_tc_efe));
    if (filterFrecuente.length) d = d.filter(r => filterFrecuente.includes(r.cliente_frecuente_si_no));
    if (filterEmpleado.length) d = d.filter(r => filterEmpleado.includes(r.empleado));
    d = d.filter(r => !r.monto || (r.monto >= filterPrecioMin && r.monto <= filterPrecioMax));
    return d;
  }, [allData, sucursal, dateMode, selectedYear, selectedMonth, exactDate, rangeStart, rangeEnd, filterSexo, filterEdad, filterPago, filterFrecuente, filterEmpleado, filterPrecioMin, filterPrecioMax]); // <--- Agrega allData a las dependencias
 
  // ── KPIs ──
  const kpis = useMemo(() => {
    const ventas = filteredData.filter(r => r.compro_si_no === "SI");
    const noVentas = filteredData.filter(r => r.compro_si_no === "NO");
    const montoTotal = _.sumBy(ventas, "monto");
    const totalClientes = _.sumBy(ventas, "num_personas");
    return {
      monto: montoTotal,
      ventas: ventas.length,
      noVentas: noVentas.length,
      clientes: totalClientes,
      ticketPromedio: ventas.length ? montoTotal / ventas.length : 0,
      productosPromedio: ventas.length ? _.meanBy(ventas, "numero_de_productos") : 0,
      conversion: filteredData.length ? (ventas.length / filteredData.length * 100) : 0,
      frecuentes: filteredData.length ? (filteredData.filter(r => r.cliente_frecuente_si_no === "SI").length / filteredData.length * 100) : 0,
    };
  }, [filteredData]);

  // ── Aggregations ──
  const getKpiValue = useCallback((rows) => {
    if (kpi === "monto") return _.sumBy(rows.filter(r => r.compro_si_no === "SI"), "monto");
    if (kpi === "ventas") return rows.filter(r => r.compro_si_no === "SI").length;
    if (kpi === "no_ventas") return rows.filter(r => r.compro_si_no === "NO").length;
    if (kpi === "clientes") return _.sumBy(rows.filter(r => r.compro_si_no === "SI"), "num_personas");
    return 0;
  }, [kpi]);

  // Ventas por hora
  const ventasHora = useMemo(() => {
    const grouped = _.groupBy(filteredData.filter(r => r.hora_num >= 9 && r.hora_num <= 21), "hora_num");
    return _.range(9, 22).map(h => ({ hora: `${String(h).padStart(2, "0")}:00`, valor: getKpiValue(grouped[h] || []) }));
  }, [filteredData, getKpiValue]);

  // Ventas por día del mes
  const ventasMes = useMemo(() => {
    const d = filteredData.filter(r => r.ano === Number(compMesYear) && r.mes_num === compMesMonth + 1);
    const grouped = _.groupBy(d, r => new Date(r.fecha).getDate());
    const daysInMonth = new Date(compMesYear, compMesMonth + 1, 0).getDate();
    return _.range(1, daysInMonth + 1).map(day => ({ dia: day, valor: getKpiValue(grouped[day] || []) }));
  }, [filteredData, compMesYear, compMesMonth, getKpiValue]);

  // Ventas por día de la semana
  const ventasSemana = useMemo(() => {
    const dias = ["Lunes", "Martes", "Miércoles", "Jueves", "Viernes", "Sábado", "Domingo"];
    const grouped = _.groupBy(filteredData, "dia_semana");
    return dias.map(d => ({ dia: d, valor: getKpiValue(grouped[d] || []) }));
  }, [filteredData, getKpiValue]);

  // Ventas anuales
  const ventasAnuales = useMemo(() => {
    const grouped = _.groupBy(filteredData, "ano");
    return Object.keys(grouped).sort().map(y => ({ anio: y, valor: getKpiValue(grouped[y]) }));
  }, [filteredData, getKpiValue]);

  // Por sucursal
  const ventasSucursal = useMemo(() => {
    const grouped = _.groupBy(filteredData, "sucursal");
    return SUCURSALES.map(s => ({ sucursal: s, valor: getKpiValue(grouped[s] || []) }));
  }, [filteredData, getKpiValue]);

  // Cliente frecuente
  const ventasFrecuente = useMemo(() => {
    const grouped = _.groupBy(filteredData, r => r.cliente_frecuente_si_no || "No especificado");
    return ["SI", "NO", "No especificado"].map(k => ({ tipo: k, valor: getKpiValue(grouped[k] || []) }));
  }, [filteredData, getKpiValue]);

  // Rango edad
  const ventasEdad = useMemo(() => {
    const orden = ["Menor", "18-25", "26-35", "36-50", "51-65", "65+", "No especificado"];
    const grouped = _.groupBy(filteredData, r => r.rango_edad || "No especificado");
    return orden.map(k => ({ rango: k, valor: getKpiValue(grouped[k] || []) }));
  }, [filteredData, getKpiValue]);

  // Tipo de pago
  const ventasPago = useMemo(() => {
    const grouped = _.groupBy(filteredData.filter(r => r.compro_si_no === "SI"), r => r.pago_tc_efe || "No especificado");
    return ["EFE", "TC", "No especificado"].map(k => ({ tipo: k === "EFE" ? "Efectivo" : k === "TC" ? "Tarjeta" : k, valor: getKpiValue(grouped[k] || []) }));
  }, [filteredData, getKpiValue]);

  // Rango precio
  const ventasPrecio = useMemo(() => {
    const rangos = [
      { label: "$1-20", min: 1, max: 20 },
      { label: "$21-100", min: 21, max: 100 },
      { label: "$101-300", min: 101, max: 300 },
      { label: "$301-1000", min: 301, max: 1000 },
      { label: "+$1000", min: 1001, max: Infinity },
    ];
    const ventas = filteredData.filter(r => r.compro_si_no === "SI" && r.monto);
    return rangos.map(rg => ({ rango: rg.label, valor: getKpiValue(ventas.filter(v => v.monto >= rg.min && v.monto <= rg.max)) }));
  }, [filteredData, getKpiValue]);

  // Top productos
  const topProductos = useMemo(() => {
    const all = [];
    filteredData.forEach(r => {
      if (!r.que_buscaba) return;
      r.que_buscaba.split(";").forEach(p => {
        const prod = p.trim().toUpperCase();
        if (prod) all.push({ ...r, producto: prod });
      });
    });
    const grouped = _.groupBy(all, "producto");
    return Object.entries(grouped)
      .map(([prod, rows]) => ({ producto: prod, valor: getKpiValue(rows) }))
      .sort((a, b) => b.valor - a.valor)
      .slice(0, 10);
  }, [filteredData, getKpiValue]);

  // Comentarios
  const topComentarios = useMemo(() => {
    const comentarios = filteredData.filter(r => r.comentarios).map(r => r.comentarios.toUpperCase().trim());
    const counted = _.countBy(comentarios);
    return Object.entries(counted).sort((a, b) => b[1] - a[1]).slice(0, 10).map(([c, n]) => ({ comentario: c, n }));
  }, [filteredData]);

  // Razones no compra
  const topRazones = useMemo(() => {
    const razones = filteredData.filter(r => r.compro_si_no === "NO" && r.razon_no_compra).map(r => r.razon_no_compra.toUpperCase().trim());
    const counted = _.countBy(razones);
    return Object.entries(counted).sort((a, b) => b[1] - a[1]).slice(0, 10).map(([r, n]) => ({ razon: r, n }));
  }, [filteredData]);

  // Empleados
  const ventasEmpleado = useMemo(() => {
    const grouped = _.groupBy(filteredData.filter(r => r.empleado), "empleado");
    return Object.entries(grouped).map(([e, rows]) => ({ empleado: e, valor: getKpiValue(rows) })).sort((a, b) => b.valor - a.valor);
  }, [filteredData, getKpiValue]);

  // Hora pico
  const horaPico = useMemo(() => {
    if (!ventasHora.length) return "--";
    const max = _.maxBy(ventasHora, "valor");
    return max ? max.hora : "--";
  }, [ventasHora]);

  // Comparison handling
  const addHoraLayer = () => {
    const name = prompt("Nombre de esta comparación:", `Comp ${horaLayers.length + 1}`);
    if (name) setHoraLayers(prev => [...prev, { name, data: ventasHora.map(h => ({ ...h })) }]);
  };
  const addMesLayer = () => {
    const name = prompt("Nombre de esta comparación:", `Comp ${mesLayers.length + 1}`);
    if (name) setMesLayers(prev => [...prev, { name, data: ventasMes.map(m => ({ ...m })) }]);
  };

  // Build comparison chart data
  const horaChartData = useMemo(() => {
    return ventasHora.map((h, i) => {
      const point = { hora: h.hora, Base: h.valor };
      horaLayers.forEach(l => { point[l.name] = l.data[i]?.valor || 0; });
      return point;
    });
  }, [ventasHora, horaLayers]);

  const mesChartData = useMemo(() => {
    return ventasMes.map((m, i) => {
      const point = { dia: m.dia, Base: m.valor };
      mesLayers.forEach(l => { point[l.name] = l.data[i]?.valor || 0; });
      return point;
    });
  }, [ventasMes, mesLayers]);

  // Table data
  const tableData = useMemo(() => {
    let d = [...allData]; // <--- AQUI CAMBIAMOS ALL_DATA por allData
    if (tableSucursal !== "Todas") d = d.filter(r => r.sucursal === tableSucursal);
    d = d.filter(r => tableCompro.includes(r.compro_si_no));
    return d;
  }, [allData, tableSucursal, tableCompro]); //

  const CHART_CATALOG = [
    { id: "ventas_hora", label: "Ventas por hora (líneas)" },
    { id: "ventas_hora_barras", label: "Ventas por hora (barras)" },
    { id: "ventas_mes", label: "Ventas del mes" },
    { id: "ventas_anual", label: "Ventas anuales" },
    { id: "ventas_semana", label: "Día de la semana" },
    { id: "cliente_frecuente", label: "Cliente frecuente" },
    { id: "rango_edad", label: "Rango de edad" },
    { id: "tipo_pago", label: "Tipo de pago" },
    { id: "rango_precio", label: "Rango de precio" },
    { id: "sucursal", label: "Comparar sucursales" },
    { id: "top_productos", label: "Top 10 productos" },
    { id: "comentarios", label: "Comentarios frecuentes" },
    { id: "motivos_no_compra", label: "Razones de no compra" },
    { id: "ventas_empleado", label: "Desempeño empleados" },
  ];

  const toggleChart = (id) => {
    setVisibleCharts(prev => prev.includes(id) ? prev.filter(c => c !== id) : [...prev, id]);
  };

  const activeFiltersCount = [filterSexo, filterEdad, filterPago, filterFrecuente, filterEmpleado].filter(f => f.length > 0).length + (filterPrecioMin > 0 || filterPrecioMax < 5000 ? 1 : 0);

  const resetFilters = () => {
    setFilterSexo([]); setFilterEdad([]); setFilterPago([]); setFilterFrecuente([]); setFilterEmpleado([]); setFilterPrecioMin(0); setFilterPrecioMax(5000);
  };

  const toggleMulti = (arr, setArr, val) => {
    setArr(prev => prev.includes(val) ? prev.filter(v => v !== val) : [...prev, val]);
  };
  if (loading) {
    return (
      <div style={{ minHeight: "100vh", background: "#0c0c1d", color: "#e2e8f0", display: "flex", justifyContent: "center", alignItems: "center", fontSize: "20px" }}>
        Cargando datos desde Railway... ⏳
      </div>
    );
  }
  // ── RENDER ──
  return (
    <div style={{ minHeight: "100vh", background: "linear-gradient(180deg, #0c0c1d 0%, #111128 50%, #0e0e20 100%)", color: "#e2e8f0", fontFamily: "'DM Sans', -apple-system, sans-serif", padding: 0 }}>
      <link href="https://fonts.googleapis.com/css2?family=DM+Sans:ital,opsz,wght@0,9..40,100..1000;1,9..40,100..1000&display=swap" rel="stylesheet" />

      {/* ── HEADER ── */}
      <header style={{ padding: "20px 28px", borderBottom: "1px solid rgba(255,255,255,0.05)", display: "flex", alignItems: "center", justifyContent: "space-between", flexWrap: "wrap", gap: 12 }}>
        <div style={{ display: "flex", alignItems: "center", gap: 14 }}>
          <div style={{ width: 38, height: 38, borderRadius: 12, background: "linear-gradient(135deg, #6366f1, #8b5cf6)", display: "flex", alignItems: "center", justifyContent: "center", fontSize: 18, fontWeight: 800 }}>V</div>
          <div>
            <h1 style={{ margin: 0, fontSize: 20, fontWeight: 800, color: "#f1f5f9", letterSpacing: "-0.5px" }}>Ventas Dashboard</h1>
            <p style={{ margin: 0, fontSize: 11, color: "#64748b" }}>Analytics en tiempo real</p>
          </div>
        </div>
        <div style={{ display: "flex", gap: 6 }}>
          <TabButton label="Resumen" active={activeTab === "resumen"} onClick={() => setActiveTab("resumen")} icon="📊" />
          <TabButton label="Datos" active={activeTab === "datos"} onClick={() => setActiveTab("datos")} icon="📋" />
        </div>
      </header>

      {activeTab === "resumen" && (
        <div style={{ padding: "20px 28px" }}>
          {/* ── FILTER BAR ── */}
          <div style={{ display: "flex", gap: 12, flexWrap: "wrap", alignItems: "flex-end", marginBottom: 20, background: "rgba(255,255,255,0.02)", borderRadius: 16, padding: "16px 20px", border: "1px solid rgba(255,255,255,0.04)" }}>
            <Select label="Sucursal" value={sucursal} onChange={setSucursal} options={["Todas", ...SUCURSALES]} />
            <Select label="Métrica" value={kpi} onChange={setKpi} options={KPI_OPTIONS} />
            <Select label="Periodo" value={dateMode} onChange={setDateMode} options={[{ value: "todo", label: "Todo" }, { value: "mes", label: "Mes" }, { value: "fecha", label: "Fecha exacta" }, { value: "rango", label: "Rango" }]} />
            {dateMode === "mes" && (
              <>
                <Select label="Año" value={selectedYear} onChange={v => setSelectedYear(Number(v))} options={[2024, 2025, 2026].map(y => ({ value: y, label: String(y) }))} />
                <Select label="Mes" value={selectedMonth} onChange={v => setSelectedMonth(Number(v))} options={MONTH_NAMES.map((m, i) => ({ value: i, label: m }))} />
              </>
            )}
            {dateMode === "fecha" && (
              <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                <label style={{ fontSize: 11, color: "#64748b", fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.5px" }}>Fecha</label>
                <input type="date" value={exactDate} onChange={e => setExactDate(e.target.value)} style={{ background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.1)", borderRadius: 8, padding: "8px 12px", color: "#e2e8f0", fontSize: 13 }} />
              </div>
            )}
            {dateMode === "rango" && (
              <>
                <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                  <label style={{ fontSize: 11, color: "#64748b", fontWeight: 600, textTransform: "uppercase" }}>Desde</label>
                  <input type="date" value={rangeStart} onChange={e => setRangeStart(e.target.value)} style={{ background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.1)", borderRadius: 8, padding: "8px 12px", color: "#e2e8f0", fontSize: 13 }} />
                </div>
                <div style={{ display: "flex", flexDirection: "column", gap: 4 }}>
                  <label style={{ fontSize: 11, color: "#64748b", fontWeight: 600, textTransform: "uppercase" }}>Hasta</label>
                  <input type="date" value={rangeEnd} onChange={e => setRangeEnd(e.target.value)} style={{ background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.1)", borderRadius: 8, padding: "8px 12px", color: "#e2e8f0", fontSize: 13 }} />
                </div>
              </>
            )}
            <button onClick={() => setShowFilters(!showFilters)} style={{ padding: "8px 16px", borderRadius: 8, border: activeFiltersCount ? "1px solid #6366f1" : "1px solid rgba(255,255,255,0.1)", background: activeFiltersCount ? "rgba(99,102,241,0.15)" : "rgba(255,255,255,0.04)", color: activeFiltersCount ? "#a5b4fc" : "#94a3b8", fontSize: 13, fontWeight: 600, cursor: "pointer", alignSelf: "flex-end", display: "flex", alignItems: "center", gap: 6 }}>
              🔍 Filtros {activeFiltersCount > 0 && <span style={{ background: "#6366f1", color: "#fff", borderRadius: "50%", width: 18, height: 18, display: "flex", alignItems: "center", justifyContent: "center", fontSize: 10, fontWeight: 800 }}>{activeFiltersCount}</span>}
            </button>
          </div>

          {/* ── ADVANCED FILTERS PANEL ── */}
          {showFilters && (
            <div style={{ background: "rgba(20,20,40,0.95)", border: "1px solid rgba(255,255,255,0.08)", borderRadius: 16, padding: "20px 24px", marginBottom: 20, animation: "fadeIn 0.2s ease" }}>
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 16 }}>
                <h3 style={{ margin: 0, fontSize: 15, fontWeight: 700, color: "#a5b4fc" }}>Filtros Avanzados</h3>
                <div style={{ display: "flex", gap: 8 }}>
                  <button onClick={resetFilters} style={{ padding: "6px 12px", borderRadius: 8, border: "none", background: "rgba(245,158,11,0.15)", color: "#fbbf24", fontSize: 12, cursor: "pointer" }}>Restablecer</button>
                  <button onClick={() => setShowFilters(false)} style={{ padding: "6px 12px", borderRadius: 8, border: "none", background: "rgba(255,255,255,0.06)", color: "#94a3b8", fontSize: 12, cursor: "pointer" }}>Cerrar</button>
                </div>
              </div>
              <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(200px, 1fr))", gap: 20 }}>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Sexo</p>
                  <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                    {["H", "M", "No especificado"].map(s => <ToggleChip key={s} label={s === "H" ? "Hombre" : s === "M" ? "Mujer" : s} active={filterSexo.includes(s)} onClick={() => toggleMulti(filterSexo, setFilterSexo, s)} />)}
                  </div>
                </div>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Rango de edad</p>
                  <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                    {[...RANGOS_EDAD, "No especificado"].map(r => <ToggleChip key={r} label={r} active={filterEdad.includes(r)} onClick={() => toggleMulti(filterEdad, setFilterEdad, r)} />)}
                  </div>
                </div>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Método de pago</p>
                  <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                    {["EFE", "TC", "No especificado"].map(p => <ToggleChip key={p} label={p === "EFE" ? "Efectivo" : p === "TC" ? "Tarjeta" : p} active={filterPago.includes(p)} onClick={() => toggleMulti(filterPago, setFilterPago, p)} />)}
                  </div>
                </div>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Cliente frecuente</p>
                  <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                    {["SI", "NO"].map(f => <ToggleChip key={f} label={f === "SI" ? "Sí" : "No"} active={filterFrecuente.includes(f)} onClick={() => toggleMulti(filterFrecuente, setFilterFrecuente, f)} />)}
                  </div>
                </div>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Empleado</p>
                  <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
                    {EMPLEADOS.map(e => <ToggleChip key={e} label={e} active={filterEmpleado.includes(e)} onClick={() => toggleMulti(filterEmpleado, setFilterEmpleado, e)} />)}
                  </div>
                </div>
                <div>
                  <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase" }}>Rango de precio</p>
                  <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
                    <input type="number" value={filterPrecioMin} onChange={e => setFilterPrecioMin(Number(e.target.value))} style={{ width: 80, background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.1)", borderRadius: 8, padding: "6px 10px", color: "#e2e8f0", fontSize: 13 }} />
                    <span style={{ color: "#64748b" }}>—</span>
                    <input type="number" value={filterPrecioMax} onChange={e => setFilterPrecioMax(Number(e.target.value))} style={{ width: 80, background: "rgba(255,255,255,0.06)", border: "1px solid rgba(255,255,255,0.1)", borderRadius: 8, padding: "6px 10px", color: "#e2e8f0", fontSize: 13 }} />
                  </div>
                </div>
              </div>
            </div>
          )}

          {/* ── KPI CARDS ── */}
          <div style={{ display: "flex", gap: 14, flexWrap: "wrap", marginBottom: 20 }}>
            <KpiCard icon="💰" label="Monto total" value={`$${kpis.monto.toLocaleString("es-MX")}`} color="#10b981" />
            <KpiCard icon="🛒" label="Ventas" value={kpis.ventas.toLocaleString()} color="#6366f1" />
            <KpiCard icon="👥" label="Clientes" value={kpis.clientes.toLocaleString()} color="#8b5cf6" />
            <KpiCard icon="❌" label="No ventas" value={kpis.noVentas.toLocaleString()} color="#ef4444" />
          </div>

          {/* ── METRICS ROW ── */}
          <div style={{ display: "flex", gap: 10, flexWrap: "wrap", marginBottom: 24 }}>
            <MetricPill icon="💰" label="Ticket promedio" value={`$${kpis.ticketPromedio.toLocaleString("es-MX", { maximumFractionDigits: 0 })}`} />
            <MetricPill icon="📦" label="Productos/venta" value={kpis.productosPromedio.toFixed(1)} />
            <MetricPill icon="✅" label="Conversión" value={`${kpis.conversion.toFixed(1)}%`} />
            <MetricPill icon="🔁" label="Frecuentes" value={`${kpis.frecuentes.toFixed(1)}%`} />
            <MetricPill icon="⏱" label="Hora pico" value={horaPico} />
          </div>

          {/* ── CHART SELECTOR ── */}
          <div style={{ marginBottom: 20 }}>
            <p style={{ fontSize: 12, color: "#64748b", fontWeight: 600, marginBottom: 8, textTransform: "uppercase", letterSpacing: "0.5px" }}>Gráficas visibles</p>
            <div style={{ display: "flex", gap: 6, flexWrap: "wrap" }}>
              {CHART_CATALOG.map(c => <ToggleChip key={c.id} label={c.label} active={visibleCharts.includes(c.id)} onClick={() => toggleChart(c.id)} />)}
            </div>
          </div>

          {/* ── CHARTS GRID ── */}
          <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(380px, 1fr))", gap: 16 }}>

            {visibleCharts.includes("ventas_hora") && (
              <ChartCard title="Ventas por hora — Líneas">
                <ComparisonManager layers={horaLayers} onAdd={addHoraLayer} onRemoveLast={() => setHoraLayers(p => p.slice(0, -1))} onClearAll={() => setHoraLayers([])} />
                <ResponsiveContainer width="100%" height={250}>
                  <LineChart data={horaChartData} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="hora" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Legend />
                    <Line type="monotone" dataKey="Base" stroke={PALETTE[0]} strokeWidth={2.5} dot={{ r: 3 }} />
                    {horaLayers.map((l, i) => <Line key={l.name} type="monotone" dataKey={l.name} stroke={PALETTE[(i + 1) % PALETTE.length]} strokeWidth={2} dot={{ r: 2 }} strokeDasharray={i % 2 ? "6 3" : undefined} />)}
                  </LineChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("ventas_hora_barras") && (
              <ChartCard title="Ventas por hora — Barras">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasHora} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="hora" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" fill={CHART_COLORS.primary} radius={[4, 4, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label} />
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("ventas_mes") && (
              <ChartCard title={`Ventas del mes — ${MONTH_NAMES[compMesMonth]} ${compMesYear}`}>
                <div style={{ display: "flex", gap: 10, marginBottom: 10, alignItems: "flex-end", flexWrap: "wrap" }}>
                  <Select label="Año" value={compMesYear} onChange={v => setCompMesYear(Number(v))} options={[2024, 2025, 2026].map(y => ({ value: y, label: String(y) }))} />
                  <Select label="Mes" value={compMesMonth} onChange={v => setCompMesMonth(Number(v))} options={MONTH_NAMES.map((m, i) => ({ value: i, label: m }))} />
                  <ComparisonManager layers={mesLayers} onAdd={addMesLayer} onRemoveLast={() => setMesLayers(p => p.slice(0, -1))} onClearAll={() => setMesLayers([])} />
                </div>
                <ResponsiveContainer width="100%" height={250}>
                  <LineChart data={mesChartData} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="dia" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Legend />
                    <Line type="monotone" dataKey="Base" stroke={PALETTE[1]} strokeWidth={2.5} dot={{ r: 2 }} />
                    {mesLayers.map((l, i) => <Line key={l.name} type="monotone" dataKey={l.name} stroke={PALETTE[(i + 2) % PALETTE.length]} strokeWidth={2} dot={{ r: 2 }} />)}
                  </LineChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("ventas_anual") && (
              <ChartCard title="Ventas anuales">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasAnuales} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="anio" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[6, 6, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasAnuales.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("ventas_semana") && (
              <ChartCard title="Ventas por día de la semana">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasSemana} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="dia" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} interval={0} angle={-25} textAnchor="end" height={50} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[6, 6, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasSemana.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("sucursal") && (
              <ChartCard title="Comparación por sucursal">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasSucursal} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="sucursal" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[6, 6, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasSucursal.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("cliente_frecuente") && (
              <ChartCard title="Cliente frecuente">
                <ResponsiveContainer width="100%" height={250}>
                  <PieChart>
                    <Pie data={ventasFrecuente} dataKey="valor" nameKey="tipo" cx="50%" cy="50%" outerRadius={90} innerRadius={50} paddingAngle={3} label={({ tipo, percent }) => `${tipo} ${(percent * 100).toFixed(0)}%`}>
                      {ventasFrecuente.map((_, i) => <Cell key={i} fill={PALETTE[i]} />)}
                    </Pie>
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                  </PieChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("rango_edad") && (
              <ChartCard title="Rango de edad">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasEdad} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="rango" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[6, 6, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasEdad.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("tipo_pago") && (
              <ChartCard title="Tipo de pago">
                <ResponsiveContainer width="100%" height={250}>
                  <PieChart>
                    <Pie data={ventasPago} dataKey="valor" nameKey="tipo" cx="50%" cy="50%" outerRadius={90} innerRadius={50} paddingAngle={3} label={({ tipo, percent }) => `${tipo} ${(percent * 100).toFixed(0)}%`}>
                      {ventasPago.map((_, i) => <Cell key={i} fill={[CHART_COLORS.primary, CHART_COLORS.secondary, "#64748b"][i]} />)}
                    </Pie>
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                  </PieChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("rango_precio") && (
              <ChartCard title="Distribución por rango de precio">
                <ResponsiveContainer width="100%" height={250}>
                  <BarChart data={ventasPrecio} margin={{ top: 10, right: 10, bottom: 0, left: 0 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="rango" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[6, 6, 0, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasPrecio.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("top_productos") && (
              <ChartCard title="Top 10 productos">
                <ResponsiveContainer width="100%" height={350}>
                  <BarChart data={topProductos} layout="vertical" margin={{ top: 0, right: 10, bottom: 0, left: 80 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis type="number" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <YAxis type="category" dataKey="producto" tick={{ fill: "#94a3b8", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} width={80} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[0, 6, 6, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {topProductos.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("comentarios") && (
              <ChartCard title="Comentarios más frecuentes">
                <ResponsiveContainer width="100%" height={350}>
                  <BarChart data={topComentarios} layout="vertical" margin={{ top: 0, right: 10, bottom: 0, left: 100 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis type="number" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis type="category" dataKey="comentario" tick={{ fill: "#94a3b8", fontSize: 10 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} width={100} />
                    <Tooltip content={<ChartTooltip />} />
                    <Bar dataKey="n" fill={CHART_COLORS.secondary} radius={[0, 6, 6, 0]} name="Frecuencia" />
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("motivos_no_compra") && (
              <ChartCard title="Razones de no compra">
                <ResponsiveContainer width="100%" height={350}>
                  <BarChart data={topRazones} layout="vertical" margin={{ top: 0, right: 10, bottom: 0, left: 110 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis type="number" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} />
                    <YAxis type="category" dataKey="razon" tick={{ fill: "#94a3b8", fontSize: 10 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} width={110} />
                    <Tooltip content={<ChartTooltip />} />
                    <Bar dataKey="n" fill={CHART_COLORS.danger} radius={[0, 6, 6, 0]} name="Frecuencia" />
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}

            {visibleCharts.includes("ventas_empleado") && (
              <ChartCard title="Desempeño por empleado">
                <ResponsiveContainer width="100%" height={280}>
                  <BarChart data={ventasEmpleado} layout="vertical" margin={{ top: 0, right: 10, bottom: 0, left: 90 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis type="number" tick={{ fill: "#64748b", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} tickFormatter={v => isMoney ? `$${(v/1000).toFixed(0)}k` : v} />
                    <YAxis type="category" dataKey="empleado" tick={{ fill: "#94a3b8", fontSize: 11 }} axisLine={{ stroke: "rgba(255,255,255,0.1)" }} width={90} />
                    <Tooltip content={<ChartTooltip isMoney={isMoney} />} />
                    <Bar dataKey="valor" radius={[0, 6, 6, 0]} name={KPI_OPTIONS.find(k => k.value === kpi)?.label}>
                      {ventasEmpleado.map((_, i) => <Cell key={i} fill={PALETTE[i % PALETTE.length]} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </ChartCard>
            )}
          </div>
        </div>
      )}

      {/* ── DATOS TAB ── */}
      {activeTab === "datos" && (
        <div style={{ padding: "20px 28px" }}>
          <div style={{ display: "flex", gap: 12, flexWrap: "wrap", marginBottom: 20, alignItems: "flex-end" }}>
            <Select label="Sucursal" value={tableSucursal} onChange={setTableSucursal} options={["Todas", ...SUCURSALES]} />
            <div>
              <p style={{ fontSize: 11, color: "#64748b", fontWeight: 600, marginBottom: 4, textTransform: "uppercase" }}>Compró</p>
              <div style={{ display: "flex", gap: 6 }}>
                {["SI", "NO"].map(v => <ToggleChip key={v} label={v} active={tableCompro.includes(v)} onClick={() => setTableCompro(p => p.includes(v) ? p.filter(x => x !== v) : [...p, v])} />)}
              </div>
            </div>
          </div>

          <div style={{ background: "rgba(20,20,40,0.8)", border: "1px solid rgba(255,255,255,0.06)", borderRadius: 16, overflow: "hidden" }}>
            <div style={{ overflowX: "auto" }}>
              <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 12 }}>
                <thead>
                  <tr style={{ borderBottom: "1px solid rgba(255,255,255,0.08)" }}>
                    {["Fecha", "Hora", "Sucursal", "Compró", "Monto", "Personas", "Productos", "Sexo", "Edad", "Pago", "Frecuente", "Empleado", "Buscaba", "Comentario"].map(h => (
                      <th key={h} style={{ padding: "12px 10px", textAlign: "left", color: "#64748b", fontWeight: 700, textTransform: "uppercase", letterSpacing: "0.5px", fontSize: 11, whiteSpace: "nowrap", position: "sticky", top: 0, background: "rgba(20,20,40,0.98)" }}>{h}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {tableData.slice(tablePage * pageSize, (tablePage + 1) * pageSize).map((r, i) => (
                    <tr key={r.id} style={{ borderBottom: "1px solid rgba(255,255,255,0.03)", background: i % 2 ? "rgba(255,255,255,0.01)" : "transparent" }}>
                      <td style={{ padding: "10px", whiteSpace: "nowrap", color: "#94a3b8" }}>{r.fecha}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.hora?.slice(0, 5)}</td>
                      <td style={{ padding: "10px", color: "#cbd5e1" }}>{r.sucursal}</td>
                      <td style={{ padding: "10px" }}><span style={{ padding: "2px 8px", borderRadius: 10, fontSize: 11, fontWeight: 600, background: r.compro_si_no === "SI" ? "rgba(16,185,129,0.15)" : "rgba(239,68,68,0.15)", color: r.compro_si_no === "SI" ? "#6ee7b7" : "#f87171" }}>{r.compro_si_no}</span></td>
                      <td style={{ padding: "10px", color: "#e2e8f0", fontWeight: 600 }}>{r.monto ? `$${r.monto.toLocaleString()}` : "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8", textAlign: "center" }}>{r.num_personas}</td>
                      <td style={{ padding: "10px", color: "#94a3b8", textAlign: "center" }}>{r.numero_de_productos || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.sexo_h_m || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.rango_edad || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.pago_tc_efe || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.cliente_frecuente_si_no}</td>
                      <td style={{ padding: "10px", color: "#94a3b8" }}>{r.empleado || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8", maxWidth: 120, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{r.que_buscaba || "—"}</td>
                      <td style={{ padding: "10px", color: "#94a3b8", maxWidth: 100, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>{r.comentarios || "—"}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", padding: "12px 16px", borderTop: "1px solid rgba(255,255,255,0.06)" }}>
              <span style={{ fontSize: 12, color: "#64748b" }}>{tableData.length.toLocaleString()} registros</span>
              <div style={{ display: "flex", gap: 6 }}>
                <button onClick={() => setTablePage(p => Math.max(0, p - 1))} disabled={tablePage === 0} style={{ padding: "6px 12px", borderRadius: 8, border: "1px solid rgba(255,255,255,0.1)", background: "rgba(255,255,255,0.04)", color: tablePage === 0 ? "#334155" : "#94a3b8", fontSize: 12, cursor: tablePage === 0 ? "default" : "pointer" }}>← Anterior</button>
                <span style={{ padding: "6px 12px", fontSize: 12, color: "#64748b" }}>Pág {tablePage + 1} de {Math.ceil(tableData.length / pageSize)}</span>
                <button onClick={() => setTablePage(p => p + 1)} disabled={(tablePage + 1) * pageSize >= tableData.length} style={{ padding: "6px 12px", borderRadius: 8, border: "1px solid rgba(255,255,255,0.1)", background: "rgba(255,255,255,0.04)", color: (tablePage + 1) * pageSize >= tableData.length ? "#334155" : "#94a3b8", fontSize: 12, cursor: (tablePage + 1) * pageSize >= tableData.length ? "default" : "pointer" }}>Siguiente →</button>
              </div>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
