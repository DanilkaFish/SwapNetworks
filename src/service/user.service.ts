import { SHA256, jwt, throwError } from "../util/utils";
import { generateJWT, generateRT } from "../util/jwt.utils";
import { PrismaClient, RefreshToken } from "@prisma/client";
import * as Models from "../model/user.model"

const prisma = new PrismaClient({});
// const SECRET_KEY = "alles"

// async function jwt_generation(rt: RefreshToken, email: string ){
//   return
// }

async function is_pw_valid(pw: string) {
    return true;
}

async function is_email_exists(email: string) {
    const user = await prisma.user.findUnique({ where: { email: email } });
    return user;
}

function makePayload(email: string, rtToken: string){
  return {
    email: email,
    iat: Date.now(),
    rt: rtToken
  }
}

function makeUserInfo(email: string, rtToken: string) {
  const payload = makePayload(email, rtToken);
    return {
        accessToken: {
            jwt: generateJWT(payload),
            payload: payload,
        },
        refreshToken: rtToken,
    };
}

export const createUser = async (email: string, password: string) => {
  const is_data_valid =
      !(await is_email_exists(email)) && (await is_pw_valid(password));
  if (!is_data_valid)
      throwError(
          "register error for : " + email + " " + password,
          400,
          "e-mail уже занят или пароль не соотвeтствует требованиям",
      );
  try {
      const user = await Models.userModel.create(email, password);
      const rt = await Models.refreshTokenModel.create(user.id);
      return makeUserInfo(email, rt.token);
  } catch (err) {
      throw err;
  }
};

export const grantAccess = async (email: string, password: string) => {
  const user = await is_email_exists(email);
  if (user) {
      if (user.password == (await SHA256(password))) {
          const rt = await Models.refreshTokenModel.updateByUserId(user.id);
          return makeUserInfo(email, rt.token);
      } else {
          throwError(
              "login error for : " + email + " " + password,
              400,
              "Пароль неверный",
          );
      }
  } else {
      throwError(
          "login error for : " + email + " " + password,
          400,
          "Пользователь с указанным e-mail не существует",
      );
  }
};

export const refreshToken = async (email: string, rtToken: string) => {
  const rt = await Models.refreshTokenModel.findByEmail(email);
  const userid = await Models.userModel.findUserIdByEmail(email);
  // const rt = 
  if (rt && userid &&(Date.now() - rt.updated_at.getTime() > 7*24*60*60*10)){
    const rt = await Models.refreshTokenModel.updateByUserId(userid);
    return makeUserInfo(email, rt.token);
  }else{
    throwError(
      "refresh token is expired : " + email + " ",
      400,
      "refresh token is expired",
    );
  }
};

// export const logout = async ()
// app.post("/", async (req, res)=> {
//   const accessToken = req.header;
//   const { header, payload, signture} = accessToken.header
//   const { email, password } = req.body;
//   const is_data_valid = (!await is_email_exists(email)) && (await is_pw_valid(password));
//   if (!is_data_valid){
//     res.status(400).json({ message: "e-mail уже занят или пароль не соответствует требованиям" });
//   } else {
//     create_user(res,email, password);
//   }
// });

// app.post("/authorization/refresh_token", async (req, res)=> {
//   const { rt } = req.body;
//   const unsignedToken = header + '.' + payload + '.';
//   const signature = HMAC-SHA256(unsignedToken, SECRET_KEY)
//   res.status(200).json( {accesToken: accessToken} )
// });
