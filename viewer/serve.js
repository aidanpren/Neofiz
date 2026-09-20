// Serves the viewer for local development: pages from viewer/, film data from runs/.
// Kept separate because runs/ is generated output and is not checked in.
//
//   ./gradlew impact           writes runs/bracket.js
//   ./gradlew smash            writes runs/smash.js
//   ./gradlew topple           writes runs/topple.js   -> /scene-viewer.html
//   ./gradlew burstFilm        writes runs/film.js     -> /burst-viewer.html
//   node viewer/serve.js       http://localhost:8731
const http = require("http"), fs = require("fs"), path = require("path");
const root = __dirname, runs = path.join(root, "..", "runs");
const types = { ".html": "text/html; charset=utf-8", ".js": "text/javascript; charset=utf-8", ".css": "text/css" };

http.createServer((req, res) => {
  let name = decodeURIComponent(req.url.split("?")[0]);
  if (name === "/") name = "/scene-viewer.html";
  const candidates = [path.join(root, name), path.join(runs, name)];
  const file = candidates.find(p => fs.existsSync(p) && fs.statSync(p).isFile());
  if (!file) {
    res.writeHead(404, { "Content-Type": "text/plain" });
    res.end(name === "/film.js"
      ? "no film yet - run ./gradlew burstFilm first"
      : name === "/bracket.js"
      ? "no film yet - run ./gradlew impact first"
      : name === "/smash.js"
      ? "no film yet - run ./gradlew smash first"
      : name === "/topple.js"
      ? "no film yet - run ./gradlew topple first"
      : "not found: " + name);
    return;
  }
  res.writeHead(200, { "Content-Type": types[path.extname(file)] || "application/octet-stream" });
  res.end(fs.readFileSync(file));
}).listen(8731, () => console.log("viewer on http://localhost:8731"));
