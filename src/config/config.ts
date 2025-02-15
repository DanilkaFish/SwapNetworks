import dotenv from "dotenv";
import fs from "fs";

const log = (message: string): void => {
    console.log(`[CONFIG]: ${message}`);
};

/**
 * Загружает файл конфигурации (.env или .env.example).
 */
export const loadConfig = (): void => {
    const envPath: string = fs.existsSync(".env") ? ".env" : ".env.example";
    dotenv.config({ path: envPath });

    log(`Configuration loaded from ${envPath}`);
};

/**
 * Извлекает и валидирует порт из переменных окружения.
 * @returns {number} Порт для сервера.
 */
export const getPort = (): number => {
    const rawPort: string | undefined = process.env.PORT;

    if (!rawPort) {
        log("PORT is not defined, using default port 3000");
        return 3000;
    }

    const parsedPort = parseInt(rawPort, 10);

    if (isNaN(parsedPort) || parsedPort <= 0 || parsedPort > 65535) {
        log(`Invalid PORT value "${rawPort}", using default port 3000`);
        return 3000;
    }

    return parsedPort;
};
