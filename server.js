// server.js (Este código corre en tu servidor, NO en el navegador)
const { Pool } = require('pg');
const express = require('express');
const cors = require('cors');

const app = express();
app.use(cors());

const pool = new Pool({
  user: 'postgres',
  host: 'switchyard.proxy.rlwy.net',
  database: 'railway',
  password: 'eNlkQFMiTIyyKoDXMQnNPPykKEXmtoep',
  port: 11813,
  ssl: { rejectUnauthorized: false } // Requerido para Railway
});

app.get('/api/datos', async (req, res) => {
  try {
    const result = await pool.query('SELECT * FROM tu_tabla_de_ventas'); 
    res.json(result.rows);
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.listen(3001, () => console.log('API lista en puerto 3001'));
