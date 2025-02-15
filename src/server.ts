import express, { Application } from "express";
import { loadConfig, getPort } from "./config/config";
import { setupRoutes } from "./route";
import { errorHandler } from "./middleware/error.middleware";
import { logger } from "./util/logger";

export const app: Application = express();
loadConfig();

const PORT: number = getPort();

app.use(express.json());
app.use(express.urlencoded({ extended: true }));
setupRoutes(app);
app.use(errorHandler);

app.listen(PORT, () => {
    logger.info(`Server running on port ${PORT}`);
});
